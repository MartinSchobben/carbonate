
##################################################################
# Steady state profiles of fast decaying organic matter 
# in aquatic sediments        
##################################################################

# required packages
library(ReacTran)  # for creating finite spatial grid and solving differential equations 
library(marelac)   # for calculating diffusion coefficients
library(emdbook)   # for generating log-transformed density distributions
library(ggplot2)   # for plotting of data
library(gridExtra) # for multigraph option


# =============================================================================
# Auxiliary functions
# =============================================================================

# ---------------------------------------
# Plotting function
# ---------------------------------------

diagen.plot <- function(depth = PL$grid$x.mid, conc,  depth.label, conc.label)
{
  plot(conc,depth,lwd=3,type="l",ylim=c(PL$L,0),
       xlim=c(min(conc), max(conc)),
       xlab="", ylab=depth.label, main = conc.label,
       axes=FALSE)
  abline(h = 0)
  axis(pos = 0, side = 2)
  axis(pos = 0, side = 3)
}

# ---------------------------
# Quality check functions
# ---------------------------

# Inventory change

IntegratedRate <- function(rate, depth = NULL)  {      # integrated rate for liquids
  if (is.null(depth))
    sum(rate * PL$grid$dx) # * PL$por.grid$mid 
  else
    sum(rate * PL$grid$dx * (PL$grid$x.mid < depth)) # * PL$por.grid$mid 
} 
IntegratedRateSolid <- function(rate, depth = NULL) {  #                     solids
  if (is.null(depth))
    sum(rate * PL$grid$dx) # * PL$svf.grid$mid 
  else
    sum(rate * PL$grid$dx * (PL$grid$x.mid < depth)) # * PL$svf.grid$mid 
} 

Sbudget <- function(output) {
  
  flux.SO4 <- output$SO4.SWI.flux - output$SO4.deep.flux
  flux.HS  <- output$HS.SWI.flux - output$HS.deep.flux
  
  flux.tot <- flux.SO4 #+ flux.HS
  
  cons.SO4 <- IntegratedRate(output$SR) + IntegratedRate(output$AOM) - IntegratedRate(output$CSO)
  cons.HS  <- - IntegratedRate(output$SR) - IntegratedRate(output$AOM) + IntegratedRate(output$CSO)    
  
  cons.tot <- cons.SO4 #+ cons.HS
  
  return(list(Total.Flux = flux.tot, SO4.Flux = flux.SO4,# HS.Flux=flux.HS,
              Total.Cons = cons.tot, SO4.Cons=cons.SO4,# HS.Cons=cons.HS, 
              Delta = flux.tot - cons.tot))
}

Cbudget <- function(output) {
  
  flux.Corg <- output$CH2O.SWI.flux - output$CH2O.deep.flux
  flux.DIC  <- output$DIC.SWI.flux - output$DIC.deep.flux
  flux.CH4  <- output$CH4.SWI.flux - output$CH4.deep.flux
  
  flux.tot <- flux.Corg + flux.DIC + flux.CH4
  
  cons.Corg  <- 2*IntegratedRateSolid(output$SR) + IntegratedRateSolid(output$AR) + IntegratedRateSolid(output$MG) 
  cons.DIC   <- - 2*IntegratedRate(output$SR) - IntegratedRate(output$AR) - 1/2*IntegratedRate(output$MG) - IntegratedRate(output$AOM) 
  cons.CH4   <- - 1/2*IntegratedRate(output$MG) + IntegratedRate(output$AOM) 
  
  cons.tot <- cons.Corg + cons.DIC + cons.CH4
  
  return(list(Total.Flux = flux.tot, Corg.Flux = flux.Corg, DIC.Flux=flux.DIC, CH4.Flux=flux.CH4,
              Total.Cons = cons.tot, Corg.Cons=cons.Corg, DIC.Cons=cons.DIC, CH4.Cons=cons.CH4, 
              Delta = flux.tot - cons.tot))
}

C12budget <- function(output) {
  
  flux.12C.CO2  <-  output$C12.CO2.SWI.flux - output$C12.CO2.deep.flux 
  flux.12C.Carb <- output$C12.Carb.SWI.flux - output$C12.Carb.deep.flux
  
  flux.tot <- flux.12C.CO2 + flux.12C.Carb
  
  cons.12C.Carb <- IntegratedRateSolid(output$Carb.12.dis) - IntegratedRateSolid(output$CO2.12.prec)
  cons.12C.CO2  <- IntegratedRate(output$CO2.12.prec) - IntegratedRate(output$Carb.12.dis) - 
    IntegratedRate(output$CO2.12.AR) - 2*IntegratedRate(output$CO2.12.SR) - 
    1/2*IntegratedRate(output$CO2.12.MG) - IntegratedRate(output$CO2.12.AOM)
  
  cons.tot <- cons.12C.Carb + cons.12C.CO2   
  
  return(list(Total.Flux = flux.tot, flux.12C.CO2 = flux.12C.CO2, flux.12C.Carb=flux.12C.Carb,
              Total.Cons = cons.tot, cons.12C.Carb=cons.12C.Carb, cons.12C.CO2=cons.12C.CO2,
              Delta = flux.tot - cons.tot))
}

C13budget <- function(output) {
  
  flux.13C.CO2  <-  output$C13.CO2.SWI.flux - output$C13.CO2.deep.flux 
  flux.13C.Carb <- output$C13.Carb.SWI.flux - output$C13.Carb.deep.flux
  
  flux.tot <- flux.13C.CO2 + flux.13C.Carb
  
  cons.13C.Carb <- IntegratedRateSolid(output$Carb.13.dis) - IntegratedRateSolid(output$CO2.13.prec)
  cons.13C.CO2  <- IntegratedRate(output$CO2.13.prec) - IntegratedRate(output$Carb.13.dis) - 
    IntegratedRate(output$CO2.13.AR) - 2*IntegratedRate(output$CO2.13.SR) - 
    1/2*IntegratedRate(output$CO2.13.MG) - IntegratedRate(output$CO2.13.AOM)
  
  cons.tot <- cons.13C.Carb + cons.13C.CO2   
  
  return(list(Total.Flux = flux.tot, flux.13C.CO2 = flux.13C.CO2, flux.13C.Carb=flux.13C.Carb,
              Total.Cons = cons.tot, cons.13C.Carb=cons.13C.Carb, cons.13C.CO2=cons.13C.CO2,
              Delta = flux.tot - cons.tot))
}

# =============================================================================
# Units used in the program: 
# =============================================================================

# Mass =  umol
# Space = cm
# Time = yr 

# =============================================================================
# Assumptions
# =============================================================================

# (1) microbial-mediated authigenic carbonate formation occurs via two prime pathways, ongoing dissolution an recrystallization with depth and as near-instantaneous authigenic carbonate formation 
# (2) the microbe's exhaled C has a static carbon isotope value specific for the metabolic pathway

PL <- list()

# =============================================================================
# Model domain and grid definition
# =============================================================================

PL$L <- 1000  # depth of sediment domain [cm]
PL$N <- 20000  # number of grid layers
PL$grid <- setup.grid.1D(x.up = 0, L = PL$L, N = PL$N)

PL$N.var <- 10 # number of state variables 13CO2, 12CO2, CH2O, O2, SO4, CH4, 12Carb, 13Carb, DIC

# =============================================================================
# Model parameters: 
# =============================================================================

# Environmental parameters

PL$S       <- 35    # salinity
PL$TC      <- 20    # temperature [deg C]
PL$P       <- 1.013 # pressure [bar]
PL$rho.s   <- 2.6   # g cm-3 solid
PL$rho.f   <- 1     # g cm-3 fluid

PL$VPDB    <- 1.1237*1E-2 # C isotope reference (13C/12C)

# =================================================================================================================================
# Physical properties of sedimentcolumn
# =================================================================================================================================

# Porosity profile 

PL$por.0  <- 0.7   # porosity
PL$por.grid <- setup.prop.1D(value = PL$por.0, grid = PL$grid)
PL$svf.grid <- setup.prop.1D(value = (1-PL$por.0), grid = PL$grid)

# Transport parameters 

PL$Db <- 5 # bioturbation intensity [cm2 yr-] Dale et al. 2016
PL$v  <- 0.2 # sedimentation rate [cm yr-]
PL$irr.0 <- 50  # Phanerozoic bioirrigation coefficient at sediment surface [yr-] Dale et al. 2016
PL$x.L <- 2 # Phanerozoic mixing depth [cm] Dale et al. 2016
PL$z.bio <- 1 # Phanerozoic attenuation depth coefficient [cm]  Dale et al. 2016

PL$Db.grid <-setup.prop.1D(func = p.exp, grid = PL$grid, y.0 = PL$Db, y.inf = 0, x.L = PL$x.L)
PL$v.grid <- setup.prop.1D(value = PL$v, grid = PL$grid)
PL$irr.grid <- setup.prop.1D(func = p.exp, grid = PL$grid, y.0 = PL$irr.0, y.inf = 0, x.att = PL$z.bio )

# =================================================================================================================================
# Reaction parameters 
# =================================================================================================================================

# The d13C isotope composition of the product (HCO3) of the following metabolic pathways; organoclastic Microbial Sulfate Reduction (MSR), Anaerobic Oxidation of Methane (AOM) and methanogenesis.

# ---------------------------------------
# Biogeochemical reactions
# ---------------------------------------

# Monod constants for metabolism limitation

PL$Ks.O2 <- 0.001 # O2 limitation for aerobic consumption of OM [umol cm-3] Faber et al., 2012 Biogeosciences
PL$Ks.SO4 <- 0.9  # SO4 limitation constant for consumption of OM by organoclastic sulfate reduction [umol cm-3] Meysman et al 2015 GCA
PL$Ks.AOM <- 1    # SO4 limitation constant for consumption of OM by methane-driven sulfate reduction [umol cm-3] Contreras et al 2013 PNAS

# first-order reaction rates for metabolism

PL$k.OM   <- 0.1    # kinetic constant for aerobic respiration [yr-1]   

PL$k.AOM  <- 1E+4   # kinetic constant for anaerobic oxidation of methane [?mol-1 cm3 yr-1] 
PL$k.CSO  <- 1E+4   # kinetic constant canonical sulfur oxidation [?mol-1 cm3 yr-1] 

# ---------------------------------------
# Isotope fractionation associated with bgc reactions
# ---------------------------------------

# Metabolism related carbon isotope fractionation factors, values taken from Irwin (1977 Nature) and Marshall (1992 Geological Magazine)

PL$AR_frac <- 1 + (-25/1000)   # aerobic organic matter remineralization
PL$SR_frac <- 1 + (-25/1000)   # organoclastic micobrial sulfate reduction
PL$AOM_frac <- 1 + (-45/1000)  # anaerobic oxidation of methane
PL$MG_frac  <- 1 + (15/1000)   # methanogenesis

# 13C abundances of metabolism produced dissolved CO2

PL$AR_ratio   <- PL$VPDB * PL$AR_frac   # aerobic organic matter remineralization
PL$SR_ratio   <- PL$VPDB * PL$SR_frac   # organoclastic micobrial sulfate reduction
PL$AOM_ratio  <- PL$VPDB * PL$AOM_frac  # anaerobic oxidation of methane
PL$MG_ratio   <- PL$VPDB * PL$MG_frac   # methanogenesis

# ---------------------------------------
# Isotope fractionation associated with carbonate precipitation 
# ---------------------------------------

# Reaction rate constants for authigenic carbonate production, Fantle & DePaolo 2007 GCA

PL$k.CCR_F <- 0.4   # calcium carbonate reactivity Fantle and Paolo 2007 GCA
PL$back.lime <-  5  # background carbon isotope value of Permian limestones relative to VPDB,  Schobben et al 2016 Chem Geol

# mineral isotope offset dissolved HCO and fixed in carbonate mineral lattice by temperature, Emrich et al. 1970 EPSL

PL$Delta_carb_bicarb <- 1.85+0.035*(PL$TC-20)
PL$back.fluid        <- PL$back.lime-PL$Delta_carb_bicarb

# Carbon isotope fractionation factor between solid (carbonate) and liquid (HCO)

PL$alpha.s <- 1 - (PL$Delta_carb_bicarb/1000)
PL$alpha.f <- 1 + (PL$Delta_carb_bicarb/1000)

# 13C abundances of dissolved HCO and fixed in carbonate mineral lattice by 20 degree Celcius

PL$Rcarb_ratio <- PL$VPDB * (PL$back.lime/1000+1)    # precipitated carbonate carbon isotope ratio Permian ocean water
PL$R_ratio     <- PL$VPDB * (1+(PL$back.fluid/1000)) # bicarbonate carbon isotope ratio Permian ocean water

# =================================================================================================================================
# Boundary conditions
# =================================================================================================================================

# ---------------------------------------
# Solids
# ---------------------------------------

# Standard flux
PL$F.up.CH2O   <- 730.5             # OC average shelf value [umol cm2- yr-1]
PL$F.up.Carb <- PL$F.up.CH2O        # carbonate flux [umol cm2- yr-1]
PL$Frac <- 0.2                      # fraction diagenetic carbonate (after Schobben et al 2016 Chem Geol)

PL$F.up.12Carb <-PL$F.up.Carb * (1/(1 + PL$Rcarb_ratio))  # [umol cm2- yr-1]
PL$F.up.13Carb <-PL$F.up.Carb * (PL$Rcarb_ratio/(1 + PL$Rcarb_ratio)) # [umol cm2- yr-1]

# ---------------------------------------
# Solutes
# ---------------------------------------

PL$C.DIC.ow    <- 4.5 # Permian DIC ocean water concentration [umol cm-3], Payne et al 2010 PNAS
PL$C.13CO2.ow  <- PL$C.DIC.ow * ((PL$C.DIC.ow * PL$R_ratio )/ (PL$C.DIC.ow  + (PL$C.DIC.ow * PL$R_ratio ))) # concentration HCO3 in overlying water [umol cm-3]
PL$C.12CO2.ow  <- PL$C.DIC.ow * (PL$C.DIC.ow/(PL$C.DIC.ow  + (PL$C.DIC.ow* PL$R_ratio ))) # concentration HCO3 in overlying water [umol cm-3]
PL$C.O2.ow     <- 0.28 # concentration O2 in overlying water [umol cm-3]
PL$C.SO4.ow    <- 4 # concentration SO4 in overlying water [umol cm-3] Permian sulfate concentration, Schobben et al 2017 Palaeo3 and Luo et al 2010 EPSL
PL$C.CH4.ow    <- 0 # concentration CH4 in overlying water [umol cm-3]
PL$C.HS.ow     <- 0 # concentration HS in overlying water [umol cm-3]

# =================================================================================================================================
# Transport parameters
# =================================================================================================================================

# ---------------------------------------
# Solutes
# ---------------------------------------

PL$tort <- 1 - 2*log(PL$por.0) # tortuosity correction

# Diffusion HCO and isotope speciation
PL$Dmol.13CO2   <- diffcoeff(S = PL$S, t = PL$TC, P = PL$P, species = "HCO3")$HCO3
PL$Dmol.13CO2   <- PL$Dmol.13CO2 * (1e4 * 3600 * 24 * 365.25) # conversion to cm2 yr-1
PL$D.13CO2      <- PL$Dmol.13CO2/PL$tort                      # tortuosity correction
PL$D.13CO2.grid <- setup.prop.1D(value = PL$D.13CO2, grid = PL$grid)

PL$Dmol.12CO2   <- diffcoeff(S = PL$S, t = PL$TC, P = PL$P, species = "HCO3")$HCO3
PL$Dmol.12CO2   <- PL$Dmol.12CO2 * (1e4 *3600 *24 * 365.25) # conversion to cm2 yr-1
PL$D.12CO2      <- PL$Dmol.12CO2/PL$tort                    # tortuosity correction
PL$D.12CO2.grid <- setup.prop.1D(value = PL$D.12CO2, grid = PL$grid)

# Diffusion O2
PL$Dmol.O2    <- diffcoeff(S = PL$S, t = PL$TC, P = PL$P, species = "O2")$O2
PL$Dmol.O2    <- PL$Dmol.O2 * (1e4 *3600 *24 * 365.25) # conversion to cm2 yr-1
PL$D.O2       <- PL$Dmol.O2/PL$tort                    # tortuosity correction
PL$D.O2.grid  <- setup.prop.1D(value = PL$D.O2, grid = PL$grid)

# Diffusion SO4
PL$Dmol.SO4   <- diffcoeff(S = PL$S, t = PL$TC, P = PL$P, species = "SO4")$SO4
PL$Dmol.SO4   <- PL$Dmol.SO4 * (1e4 *3600 *24 * 365.25) # conversion to cm2 yr-1
PL$D.SO4      <- PL$Dmol.SO4/PL$tort                    # tortuosity correction
PL$D.SO4.grid <- setup.prop.1D(value = PL$D.SO4, grid = PL$grid)

# Diffusion CH4
PL$Dmol.CH4   <- diffcoeff(S = PL$S, t = PL$TC, P = PL$P, species = "CH4")$CH4
PL$Dmol.CH4   <- PL$Dmol.CH4 * (1e4 *3600 *24 * 365.25) # conversion to cm2 yr-1
PL$D.CH4      <- PL$Dmol.CH4/PL$tort                    # tortuosity correction
PL$D.CH4.grid <- setup.prop.1D(value = PL$D.CH4, grid = PL$grid)

# Diffusion DIC
PL$Dmol.DIC   <- diffcoeff(S = PL$S, t = PL$TC, P = PL$P, species = "HCO3")$HCO3
PL$Dmol.DIC   <- PL$Dmol.DIC * (1e4 *3600 *24 * 365.25) # conversion to cm2 yr-1
PL$D.DIC      <- PL$Dmol.DIC/PL$tort                    # tortuosity correction
PL$D.DIC.grid <- setup.prop.1D(value = PL$D.DIC, grid = PL$grid)

# Diffusion CH4
PL$Dmol.HS    <- diffcoeff(S = PL$S, t = PL$TC, P = PL$P, species = "H2S")$H2S
PL$Dmol.HS    <- PL$Dmol.HS * (1e4 *3600 *24 * 365.25) # conversion to cm2 yr-1
PL$D.HS       <- PL$Dmol.HS/PL$tort                    # tortuosity correction
PL$D.HS.grid  <- setup.prop.1D(value = PL$D.HS, grid = PL$grid)



# =============================================================================
# Model formulation
# =============================================================================

multi.model <- function (t,state,parameters) 
{
  with(as.list(c(parameters)),{
    
    # -----------------------------------
    # Initialisation of state variables 
    # -----------------------------------
    
    C.12CO2   <- state[1:N]
    C.13CO2   <- state[(N+1):(2*N)]
    C.O2      <- state[(2*N+1):(3*N)]
    C.SO4     <- state[(3*N+1):(4*N)]
    C.CH4     <- state[(4*N+1):(5*N)]
    C.CH2O    <- state[(5*N+1):(6*N)]
    C.12Carb  <- state[(6*N+1):(7*N)]
    C.13Carb  <- state[(7*N+1):(8*N)]
    C.DIC     <- state[(8*N+1):(9*N)]
    C.HS      <- state[(9*N+1):(10*N)]
    
    # -----------------------------------
    # Transport terms
    # -----------------------------------
    
    # Solutes
    
    Tran.12CO2 <- tran.1D(C=C.12CO2 , C.up=C.12CO2.ow , D=D.12CO2.grid, v=v.grid , VF=por.grid , dx=grid,full.check=TRUE)$dC
    Tran.13CO2 <- tran.1D(C=C.13CO2 , C.up=C.13CO2.ow , D=D.13CO2.grid, v=v.grid , VF=por.grid , dx=grid,full.check=TRUE)$dC
    Tran.O2    <- tran.1D(C=C.O2 , C.up=C.O2.ow , D=D.O2.grid , VF=por.grid, v=v.grid , dx=grid,full.check=TRUE)$dC
    Tran.SO4   <- tran.1D(C=C.SO4 , C.up=C.SO4.ow , D=D.SO4.grid, v=v.grid , VF=por.grid , dx=grid,full.check=TRUE)$dC
    Tran.CH4   <- tran.1D(C=C.CH4 , C.up=C.CH4.ow , D=D.CH4.grid, v=v.grid , VF=por.grid , dx=grid,full.check=TRUE)$dC
    Tran.DIC   <- tran.1D(C=C.DIC , C.up=C.DIC.ow , D=D.DIC.grid, v=v.grid , VF=por.grid , dx=grid,full.check=TRUE)$dC
    Tran.HS    <- tran.1D(C=C.HS , C.up=C.HS.ow , D=D.HS.grid, v=v.grid , VF=por.grid , dx=grid,full.check=TRUE)$dC
    
    # Solids
    
    Tran.CH2O   <- tran.1D(C=C.CH2O, flux.up = F.up.CH2O , v=v.grid , VF = svf.grid , dx=grid,full.check=TRUE)$dC
    Tran.12Carb <- tran.1D(C=C.12Carb, flux.up = F.up.12Carb*Frac , v=v.grid , VF = svf.grid , dx=grid,full.check=TRUE)$dC
    Tran.13Carb <- tran.1D(C=C.13Carb, flux.up = F.up.13Carb*Frac , v=v.grid , VF = svf.grid , dx=grid,full.check=TRUE)$dC
    
    
    # non-local transport term for solutes
    
    Irr.12CO2  <- irr.grid$mid*(C.12CO2.ow - C.12CO2)
    Irr.13CO2  <- irr.grid$mid*(C.13CO2.ow - C.13CO2)
    Irr.O2     <- irr.grid$mid*(C.O2.ow    - C.O2)
    Irr.SO4    <- irr.grid$mid*(C.SO4.ow   - C.SO4)
    Irr.CH4    <- irr.grid$mid*(C.CH4.ow   - C.CH4)
    Irr.DIC    <- irr.grid$mid*(C.DIC.ow   - C.DIC)
    Irr.HS     <- irr.grid$mid*(C.HS.ow    - C.HS)
    
    
    # -----------------------------------
    # Reaction terms
    # -----------------------------------
    
    # Organic matter mineralisation (OMM)
    
    Cmin <- svf.grid$mid*k.OM*C.CH2O
    
    O2.lim  <- (C.O2/(Ks.O2+C.O2))
    O2.inh  <- (Ks.O2/(C.O2+Ks.O2))
    SO4.lim <- (C.SO4/(Ks.SO4+C.SO4))
    SO4.inh <- (Ks.SO4/(C.SO4+Ks.SO4))
    
    a.f <- O2.lim/(O2.lim + SO4.lim*O2.inh + SO4.inh*O2.inh)*(C.O2>0)
    s.f <- (SO4.lim*O2.inh)/(O2.lim + SO4.lim*O2.inh + SO4.inh*O2.inh)*(C.SO4>0)
    m.f <- (SO4.inh*O2.inh)/(O2.lim + SO4.lim*O2.inh + SO4.inh*O2.inh)
    
    AR <- a.f*Cmin
    SR <- s.f*Cmin
    MG <- m.f*Cmin
    
    # Isotope fractionation associated with OMM
    
    CO2.12.AR <- AR*(1/(1  + AR_ratio))
    CO2.13.AR <- AR*(AR_ratio/(1  + AR_ratio)) 
    
    CO2.12.SR <- SR*(1/(1 + SR_ratio))
    CO2.13.SR <- SR*(SR_ratio/(1 + SR_ratio))   
    
    CO2.12.MG <- MG*(1/(1  + MG_ratio))   
    CO2.13.MG <- MG*(MG_ratio/(1  + MG_ratio))   
    
    # Re-oxidation of reduced species
    
    AOM <- por.grid$mid*k.AOM*C.SO4*C.CH4*SO4.lim*(C.CH4>0)*(C.SO4>0)
    CSO <- por.grid$mid*k.CSO*C.HS*C.O2*SO4.lim*(C.HS>0)*(C.O2>0)
    
    # Isotope fractionation associated with AOM
    
    CO2.12.AOM <- AOM*(1/(1  + AOM_ratio))   
    CO2.13.AOM <- AOM*(AOM_ratio/(1  + AOM_ratio))  
    
    # -----------------------------------
    # Precipitation – re-precipitation isotope fractionation
    # -----------------------------------
    
    Carb.dis <- (k.CCR_F/Frac)*exp(-(((PL$grid$x.mid/PL$v)/1E6)/(0.876*Frac)))
    CO2.prec <- (k.CCR_F/Frac)*exp(-(((PL$grid$x.mid/PL$v)/1E6)/(0.876*Frac)))
    
    Carb.12.dis <- Carb.dis*(C.12Carb/(C.12Carb+alpha.s*C.13Carb)) * O2.inh
    Carb.13.dis <- Carb.dis*(alpha.s*C.13Carb/(C.12Carb+alpha.s*C.13Carb)) *O2.inh
    
    CO2.12.prec <- CO2.prec*(C.12CO2/(C.12CO2+alpha.f*C.13CO2))* O2.inh
    CO2.13.prec <- CO2.prec*(C.13CO2*alpha.f/(C.12CO2+alpha.f*C.13CO2))* O2.inh
    
    # -----------------------------------
    # Combine reaction terms 
    # -----------------------------------
    
    Reac.O2    <- (- AR - 2*CSO)/por.grid$mid
    Reac.SO4   <- (- SR - AOM + CSO)/por.grid$mid
    Reac.CH4   <- (+ 1/2*MG - AOM)/por.grid$mid
    Reac.DIC   <- (AR + 2*SR + 1/2*MG + AOM)/por.grid$mid
    Reac.HS    <- (+ SR + AOM - CSO)/por.grid$mid
    
    Reac.12CO2 <- (CO2.12.AR + 2*CO2.12.SR + 1/2*CO2.12.MG + CO2.12.AOM - CO2.12.prec + Carb.12.dis)/por.grid$mid
    Reac.13CO2 <- (CO2.13.AR + 2*CO2.13.SR + 1/2*CO2.13.MG + CO2.13.AOM - CO2.13.prec + Carb.13.dis)/por.grid$mid
    
    Reac.CH2O   <- (- AR - 2*SR - MG)/svf.grid$mid
    
    
    Reac.12Carb <- (- Carb.12.dis + CO2.12.prec)/svf.grid$mid
    Reac.13Carb <- (- Carb.13.dis + CO2.13.prec)/svf.grid$mid
    
    # -----------------------------------
    # Partial differential equations
    # -----------------------------------
    
    ddt.O2    <- Tran.O2    + Reac.O2     + Irr.O2
    ddt.SO4   <- Tran.SO4   + Reac.SO4    + Irr.SO4
    ddt.CH4   <- Tran.CH4   + Reac.CH4    + Irr.CH4
    ddt.HS    <- Tran.HS    + Reac.HS     + Irr.HS
    ddt.DIC   <- Tran.DIC   + Reac.DIC    + Irr.DIC
    ddt.12CO2 <- Tran.12CO2 + Reac.12CO2  + Irr.12CO2
    ddt.13CO2 <- Tran.13CO2 + Reac.13CO2  + Irr.13CO2
    
    
    ddt.CH2O   <- Tran.CH2O + Reac.CH2O 
    ddt.12Carb <- Tran.12Carb + Reac.12Carb
    ddt.13Carb <- Tran.13Carb + Reac.13Carb
    
    # Assemble the total rate of change
    SO4.SWI.flux   <- tran.1D(C=C.SO4,C.up=C.SO4.ow,D=D.SO4.grid,v=v.grid,VF=por.grid,dx=grid)$flux.up 
    SO4.deep.flux  <- tran.1D(C=C.SO4,C.up=C.SO4.ow,D=D.SO4.grid,v=v.grid,VF=por.grid,dx=grid)$flux.down
    DIC.SWI.flux   <- tran.1D(C=C.DIC,C.up=C.DIC.ow,D=D.DIC.grid,v=v.grid,VF=por.grid,dx=grid)$flux.up
    DIC.deep.flux  <- tran.1D(C=C.DIC,C.up=C.DIC.ow,D=D.DIC.grid,v=v.grid,VF=por.grid,dx=grid)$flux.down
    HS.SWI.flux    <- tran.1D(C=C.HS,C.up=C.HS.ow,D=D.HS.grid,v=v.grid,VF=por.grid,dx=grid)$flux.up         
    HS.deep.flux   <- tran.1D(C=C.HS,C.up=C.HS.ow,D=D.HS.grid,v=v.grid,VF=por.grid,dx=grid)$flux.down    
    O2.SWI.flux    <- tran.1D(C=C.O2,C.up=C.O2.ow,D=D.O2.grid,VF=por.grid,v=v.grid,dx=grid)$flux.up         
    O2.deep.flux   <- tran.1D(C=C.O2,C.up=C.O2.ow,D=D.O2.grid,VF=por.grid,v=v.grid,dx=grid)$flux.down    
    CH4.SWI.flux   <- tran.1D(C=C.CH4,C.up=C.CH4.ow,D=D.CH4.grid,VF=por.grid,v=v.grid,dx=grid)$flux.up         
    CH4.deep.flux  <- tran.1D(C=C.CH4,C.up=C.CH4.ow,D=D.CH4.grid,VF=por.grid,v=v.grid,dx=grid)$flux.down    
    
    C12.CO2.SWI.flux  <- tran.1D(C=C.12CO2,C.up=C.12CO2.ow,D=D.12CO2.grid,v=v.grid,VF=por.grid,dx=grid)$flux.up
    C12.CO2.deep.flux <- tran.1D(C=C.12CO2,C.up=C.12CO2.ow,D=D.12CO2.grid,v=v.grid,VF=por.grid,dx=grid)$flux.down
    C13.CO2.SWI.flux  <- tran.1D(C=C.13CO2,C.up=C.13CO2.ow,D=D.13CO2.grid,v=v.grid,VF=por.grid,dx=grid)$flux.up
    C13.CO2.deep.flux <- tran.1D(C=C.13CO2,C.up=C.13CO2.ow,D=D.13CO2.grid,v=v.grid,VF=por.grid,dx=grid)$flux.down
    
    CH2O.SWI.flux    <- tran.1D(C=C.CH2O,flux.up=F.up.CH2O,v=v.grid,VF=svf.grid,dx=grid)$flux.up        
    CH2O.deep.flux   <- tran.1D(C=C.CH2O,flux.up=F.up.CH2O,v=v.grid,VF=svf.grid,dx=grid)$flux.down      
    
    C12.Carb.SWI.flux  <- tran.1D(C=C.12Carb,flux.up=F.up.12Carb,v=v.grid,VF=svf.grid,dx=grid)$flux.up 
    C12.Carb.deep.flux <- tran.1D(C=C.12Carb,flux.up=F.up.12Carb,v=v.grid,VF=svf.grid,dx=grid)$flux.down 
    C13.Carb.SWI.flux  <- tran.1D(C=C.13Carb,flux.up=F.up.13Carb,v=v.grid,VF=svf.grid,dx=grid)$flux.up 
    C13.Carb.deep.flux <- tran.1D(C=C.13Carb,flux.up=F.up.13Carb,v=v.grid,VF=svf.grid,dx=grid)$flux.down 
    
    return(list(c(ddt.12CO2 , ddt.13CO2 , ddt.O2 , ddt.SO4 , ddt.CH4, ddt.CH2O, ddt.12Carb, ddt.13Carb, ddt.DIC, ddt.HS), 
                
                # Reactions
                
                AR=AR,  SR=SR, MG=MG, AOM=AOM, CSO=CSO,
                
                Carb.12.dis=Carb.12.dis, Carb.13.dis = Carb.13.dis,
                CO2.12.prec=CO2.12.prec, CO2.13.prec=CO2.13.prec,
                CO2.12.AR=CO2.12.AR,CO2.13.AR=CO2.13.AR,
                CO2.12.SR=CO2.12.SR,CO2.13.SR=CO2.13.SR,
                CO2.12.MG=CO2.12.MG,CO2.13.MG=CO2.13.MG,
                CO2.12.AOM=CO2.12.AOM,CO2.13.AOM=CO2.13.AOM,
                
                # Fluxes
                
                SO4.SWI.flux=SO4.SWI.flux,SO4.deep.flux=SO4.deep.flux,
                DIC.SWI.flux=DIC.SWI.flux,DIC.deep.flux=DIC.deep.flux,
                HS.SWI.flux=HS.SWI.flux,HS.deep.flux=HS.deep.flux,
                O2.SWI.flux=O2.SWI.flux,O2.deep.flux=O2.deep.flux,
                CH4.SWI.flux=CH4.SWI.flux,CH4.deep.flux=CH4.deep.flux,
                C12.CO2.SWI.flux=C12.CO2.SWI.flux,C12.CO2.deep.flux=C12.CO2.deep.flux,
                C13.CO2.SWI.flux=C13.CO2.SWI.flux,C13.CO2.deep.flux=C13.CO2.deep.flux,
                
                CH2O.SWI.flux=CH2O.SWI.flux,CH2O.deep.flux=CH2O.deep.flux,
                C12.Carb.SWI.flux=C12.Carb.SWI.flux,C12.Carb.deep.flux=C12.Carb.deep.flux,
                C13.Carb.SWI.flux=C13.Carb.SWI.flux,C13.Carb.deep.flux=C13.Carb.deep.flux
                
    ))
    
  })}  

# =============================================================================
# Model solution     
# =============================================================================

# ---------------------------------------
# Model test
# ---------------------------------------

# Initial conditions

C.12CO2.in   <- rep(PL$C.12CO2.ow,length.out = PL$N)
C.13CO2.in   <- rep(PL$C.13CO2.ow,length.out = PL$N)
C.O2.in      <- rep(PL$C.O2.ow,length.out = PL$N)
C.SO4.in     <- rep(PL$C.SO4.ow,length.out = PL$N)
C.CH4.in     <- rep(PL$C.CH4.ow,length.out = PL$N)
C.CH2O.in    <- rep(1, length.out = PL$N)
C.12Carb.in  <- rep(1, length.out = PL$N)
C.13Carb.in  <- rep(1, length.out = PL$N)
C.DIC.in     <- rep(PL$C.DIC.ow, length.out = PL$N)
C.HS.in      <- rep(PL$C.HS.ow,  length.out = PL$N)

# Initialization state variables vector 

state <- c(C.12CO2.in, C.13CO2.in , C.O2.in, C.SO4.in  , C.CH4.in , C.CH2O.in,  C.12Carb.in , C.13Carb.in, C.DIC.in, C.HS.in)





#================================================================================================================================
# Plotting of model solutions, diagenetic depth profiles
#================================================================================================================================
# Plotting of redox zonation and accompanying biochemical reactions as depth profiles, in terms of prime pore water solutes and carbon isotope composition of both DIC and the solid carbonate fraction 
#----------------------------------------------------------------------

theme=theme_set(theme_classic()) # global environment


# Redox zonation with average OC accumulation, 500 [umol cm2- yr-1]

PL$F.up.CH2O <- 500

# Steady state: numerical solution

output <- steady.1D(y = state, func = multi.model, parms = PL,
                    nspec = PL$N.var, positive = TRUE )

# Transformation of model outcome 

Sens.Pore<-as.data.frame(matrix(NA, nrow=PL$N, ncol=9))
colnames(Sens.Pore)<-c("depth", "d13CDIC","DIC", "O2", "SO4", "CH4", "OC", "d13Ccarb", "HS" )

Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),1] <- PL$grid$x.mid
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),2] <- (((output$y[(1*PL$N+1):(2*PL$N)]/output$y[1:PL$N])/(PL$VPDB))-1)*1000 # del13C bicarbonate
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),3] <- output$y[(8*PL$N+1):(9*PL$N)] #DIC
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),4] <- output$y[(2*PL$N+1):(3*PL$N)] #O2
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),5] <- output$y[(3*PL$N+1):(4*PL$N)] #SO4
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),6] <- output$y[(4*PL$N+1):(5*PL$N)] #CH4
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),7] <- output$y[(5*PL$N+1):(6*PL$N)] #OC
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),8] <- (((((output$y[(7*PL$N+1):(8*PL$N)]/output$y[(6*PL$N+1):(7*PL$N)])/(PL$VPDB))-1)*1000)*PL$Frac) + ((PL$back.lime)*(1-PL$Frac))  # del13C carbonate
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),9] <- output$y[(9*PL$N+1):(10*PL$N)] #HS

# Average OC flux diagenetic depth profile plot

A<-ggplot(Sens.Pore, aes(y=depth, x=O2))+ 
  scale_y_reverse(limits=c(10,0))+
  geom_rect(aes(xmax=-0.6, xmin=-0.1,  ymax=0, ymin=2), fill="grey")+
  geom_rect(aes(xmax=-0.6, xmin=-0.1,  ymax=2, ymin=10), fill="aquamarine1")+
  geom_point(shape=1, color="#006d2c")+
  geom_point(aes(x=SO4), shape=1, color="#b2e2e2")+ 
  geom_point(aes(x=HS), shape=1, color="#66c2a4")+ 
  geom_point(aes(x=CH4), shape=1, color="#2ca25f")+ 
  
  ylab("depth (cm)")+
  xlab(expression(atop(paste(mu*mol~cm^-3))))+
  annotate("text", label=c("microbial-sulfate reduction"), y=5,  x=-0.35, size=3.5, angle=90)+
  annotate("text", label="O2", y=0,  x=0.75, size=3.5)+
  annotate("text", label="CH4", y=2,  x=0.6 , size=3.5)+
  annotate("text", label="HS", y=5,  x=1, size=3.5)+
  annotate("text", label="SO4", y=5,  x=2.5, size=3.5)+
  
  ggtitle("(a)"~~"average OC")+
  theme(
    axis.line.x = element_line(color = "black"), 
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust=(0)))

# d13C profile with average OC accumulation, 500 [umol cm2- yr-1] high OC flux d13C profile plot

B<-ggplot(Sens.Pore, aes(y=depth, x=d13CDIC))+ 
  scale_y_reverse(limits=c(100,0))+
  geom_rect(aes(xmax=-16, xmin=-14,  ymax=0, ymin=2), fill="grey")+
  geom_rect(aes(xmax=-16, xmin=-14,  ymax=2, ymin=10), fill="aquamarine1")+
  geom_point(shape=1, color="gray48")+
  geom_point(aes(x=d13Ccarb), shape=1, color="cyan3")+  
  xlim(-16,16)+
  ylab("depth (cm)")+
  xlab(expression(atop(paste(delta^13*C~"(\211 VPDB)"))))+
  annotate("text", label=c("microbial-sulfate reduction >"), y=60,  x=-15, size=3.5, angle=90)+
  annotate("text", label="DIC", y=15,  x=-9.5, size=3.5)+
  annotate("text", label="carbonate", y=30,  x=6.5 , size=3.5, angle=90)+
  
  
  ggtitle("(b)"~~"average OC")+
  theme(
    axis.line.x = element_line(color = "black"), 
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust=(0)))

# Redox zonation with average OC accumulation, 1200 [umol cm2- yr-1]

PL$F.up.CH2O <- 1200

# Steady state: numerical solution

output <- steady.1D(y = state, func = multi.model, parms = PL,
                    nspec = PL$N.var, positive = TRUE )

# Transformation of model outcome 

Sens.Pore<-as.data.frame(matrix(NA, nrow=PL$N, ncol=9))
colnames(Sens.Pore)<-c("depth", "d13CDIC","DIC", "O2", "SO4", "CH4", "OC", "d13Ccarb", "HS" )

Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),1] <- PL$grid$x.mid
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),2] <- (((output$y[(1*PL$N+1):(2*PL$N)]/output$y[1:PL$N])/(PL$VPDB))-1)*1000 # del13C bicarbonate
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),3] <- output$y[(8*PL$N+1):(9*PL$N)] #DIC
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),4] <- output$y[(2*PL$N+1):(3*PL$N)] #O2
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),5] <- output$y[(3*PL$N+1):(4*PL$N)] #SO4
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),6] <- output$y[(4*PL$N+1):(5*PL$N)] #CH4
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),7] <- output$y[(5*PL$N+1):(6*PL$N)] #OC
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),8] <- (((((output$y[(7*PL$N+1):(8*PL$N)]/output$y[(6*PL$N+1):(7*PL$N)])/(PL$VPDB))-1)*1000)*PL$Frac) + ((PL$back.lime)*(1-PL$Frac))  # del13C carbonate
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),9] <- output$y[(9*PL$N+1):(10*PL$N)] #HS  

# High OC flux diagenetic depth profile plot

C<-ggplot(Sens.Pore, aes(y=depth, x=O2))+ 
  scale_y_reverse(limits=c(10,0))+
  geom_rect(aes(xmax=-0.6, xmin=-0.1,  ymax=0, ymin=0.5), fill="aquamarine1")+
  geom_rect(aes(xmax=-0.6, xmin=-0.1,  ymax=0.5, ymin=2.5), fill="aquamarine3")+
  geom_rect(aes(xmax=-0.6, xmin=-0.1,  ymax=2.5, ymin=10), fill="lightgreen")+
  geom_point(shape=1, color="#006d2c")+
  geom_point(aes(x=SO4), shape=1, color="#b2e2e2")+ 
  geom_point(aes(x=HS), shape=1, color="#66c2a4")+ 
  geom_point(aes(x=CH4), shape=1, color="#2ca25f")+ 
  
  ylab("depth (cm)")+
  xlab(expression(atop(paste(mu*mol~cm^-3))))+
  
  annotate("text", label=c("AOM"), y=1.5,  x=-0.35, size=3.5, angle=90)+
  annotate("text", label=c("organoclastic methanogenesis"), y=6,  x=-0.35, size=3.5, angle=90)+
  annotate("text", label="CH4", y=7.5,  x=1.2 , size=3.5)+
  annotate("text", label="HS", y=5,  x=3 , size=3.5)+
  annotate("text", label="SO4", y=1.8,  x=1, size=3.5)+
  
  ggtitle("(c)"~~"high OC")+
  theme(
    axis.line.x = element_line(color = "black"), 
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust=(0)))



# d13C profile with average OC accumulation, 7000 [umol cm2- yr-1]

PL$F.up.CH2O <- 7000

# Steady state: numerical solution

output <- steady.1D(y = state, func = multi.model, parms = PL,
                    nspec = PL$N.var, positive = TRUE )

# Transformation of model outcome 

Sens.Pore<-as.data.frame(matrix(NA, nrow=PL$N, ncol=9))
colnames(Sens.Pore)<-c("depth", "d13CDIC","DIC", "O2", "SO4", "CH4", "OC", "d13Ccarb", "HS" )

Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),1] <- PL$grid$x.mid
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),2] <- (((output$y[(1*PL$N+1):(2*PL$N)]/output$y[1:PL$N])/(PL$VPDB))-1)*1000 # del13C bicarbonate
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),3] <- output$y[(8*PL$N+1):(9*PL$N)] #DIC
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),4] <- output$y[(2*PL$N+1):(3*PL$N)] #O2
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),5] <- output$y[(3*PL$N+1):(4*PL$N)] #SO4
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),6] <- output$y[(4*PL$N+1):(5*PL$N)] #CH4
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),7] <- output$y[(5*PL$N+1):(6*PL$N)] #OC
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),8] <- (((((output$y[(7*PL$N+1):(8*PL$N)]/output$y[(6*PL$N+1):(7*PL$N)])/(PL$VPDB))-1)*1000)*PL$Frac) + ((PL$back.lime)*(1-PL$Frac))  # del13C carbonate
Sens.Pore[c((((PL$N)-PL$N)+1):(PL$N)),9] <- output$y[(9*PL$N+1):(10*PL$N)] #HS  

# High OC flux d13C profile plot

D<-ggplot(Sens.Pore, aes(y=depth, x=d13CDIC))+ 
  scale_y_reverse(limits=c(100,0))+
  geom_rect(aes(xmax=-16, xmin=-14,  ymax=0, ymin=0.5), fill="aquamarine1")+
  geom_rect(aes(xmax=-16, xmin=-14,  ymax=0.5, ymin=2.5), fill="aquamarine3")+
  geom_rect(aes(xmax=-16, xmin=-14,  ymax=2.5, ymin=10), fill="lightgreen")+
  geom_point(shape=1, color="gray48")+
  geom_point(aes(x=d13Ccarb), shape=1, color="cyan3")+  
  
  ylab("depth (cm)")+
  xlim(-16,16)+
  xlab(expression(atop(paste(delta^13*C~"(\211 VPDB)"))))+
  annotate("text", label=c("organoclastic methanogenesis >"), y=60,  x=-15, size=3.5, angle=90)+
  annotate("text", label="DIC", y=15,  x=12, size=3.5)+
  annotate("text", label="carbonate", y=30,  x=3.2 , size=3.5, angle=90)+
  
  
  ggtitle("(d)"~~"high OC" )+
  theme(
    axis.line.x = element_line(color = "black"), 
    axis.line.y = element_line(color = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust=(0)))



grid.arrange(A,B,C,D, ncol=4)

depth<-arrangeGrob(A,B,C, D, ncol=4,nrow=1)

ggsave("Figure4.tiff",depth, height= 12, width=24, units="cm" ) # saves produced graphs


# ============================================================================================================================
# Sensitivity tests sedimentary and oceanic parameters
# ============================================================================================================================


# Modern OC accumulation for model calibration, data from Muller and Suess (1979) Deep Sea Research Part A 

org.acc    <- c(0.16*1E4/12 ,0.055*1E4/12 , 0.327*1E4/12 , 1.20*1E4/12 , 0.076*1E4/12 , 0.066*1E4/12 ,0.332*1E4/12 , 1.21*1E4/12 , 1.53*1E4/12 , 1.08*1E4/12 , 0.058*1E4/12 , 0.0020*1E4/12 , 0.0068*1E4/12, 0.0051*1E4/12 , 0.0042*1E4/12 ,0.0036*1E4/12 , 0.0052*1E4/12 , 0.0037*1E4/12 , 39.7*1E4/12 , 37.8*1E4/12 ,0.436*1E4/12 ) # real data org C accumulation 

org.acc2    <- seq(min(org.acc),max(org.acc) ,100) # observed natural range of OC flux 


S.P<-list()
S.P$Frac      <- c(0.1, 0.2,0.3,0.4)         # sensitivity parameter authigenic carbonate fraction  
S.P$O2.ow     <- c(0.001, 0.01, 0.1, 0.28)   # sensitivity parameter oxygen level overlying watercolumn  
S.P$SO4.ow    <- c(1, 4, 10, 28)             # sensitivity parameter sulfate level overlying watercolumn
S.P$DIC.ow    <- seq(2,7, length.out=4)      # sensitivity parameter DIC overlying watercolumn, Phanerozoic range from Ridgwell, 2005 
S.P$sed.rate  <- c( 0.2 ,0.3, 0.4, 0.5)      # sensitivity parameter sedimentation rate 
Out           <- list(NA,NA,NA,NA,NA)        # empty list to collect results of sensitivity tests

#-------------------------------------------------------------------------------------------------------------------
# Sensitivity test run sedimentary and oceanic parameters
#-------------------------------------------------------------------------------------------------------------------

system.time({
  
  for(j in c(1:length(S.P))){
    
    
    var<-S.P[[j]]
    
    # summary data matrix
    Sens.sum<-matrix(NA, ncol= 6, nrow= length(org.acc2)*length(var))
    colnames(Sens.sum)<-c("d13Cdia", "D13Cdia", "d13Cauth", "D13Cauth","var.sense", "FOC" )
    Sens.sum<-as.data.frame(Sens.sum)
    
    for(k in c(1:length(var))){
      
      for(t in c(1:length(org.acc2))){
        
        PL$F.up.CH2O<- org.acc2[t] # the OC accumulation flux
        
        #-------------------------------------------------------------------------------------------------------------------
        # Sensitivity of the model to the following changing variables with fraction of authigenic carbonate, dissolved oxygen concentration of the overlying watercolumn, dissolved sulfate concentration and DIC of the overlying watercolumn and sedimentation rate. The chosen values are considered  to be representative of a latest Permian carbonate depositional setting.
        #-------------------------------------------------------------------------------------------------------------------   
        
        if(names(S.P)[j]=="Frac"){PL$Frac<-var[k] } else { PL$Frac<- 0.2 }         # fraction of authigenic carbonate
        if(names(S.P)[j]=="O2.ow"){PL$C.O2.ow <-var[k]} else { PL$C.O2.ow <- 0.28} # oxygen level
        if(names(S.P)[j]=="SO4.ow"){PL$C.SO4.ow <- var[k]} else{ PL$C.SO4.ow <- 4} # sulfate level
        if(names(S.P)[j]=="DIC.ow"){PL$C.DIC.ow <- var[k]
        
        PL$C.13CO2.ow  <- PL$C.DIC.ow * ((PL$C.DIC.ow * PL$R_ratio )/ (PL$C.DIC.ow  + (PL$C.DIC.ow * PL$R_ratio ))) # concentration HCO3 in overlying water [umol cm-3]
        PL$C.12CO2.ow  <- PL$C.DIC.ow * (PL$C.DIC.ow/(PL$C.DIC.ow  + (PL$C.DIC.ow* PL$R_ratio ))) # concentration HCO3 in overlying water [umol cm-3]
        }else{ PL$C.DIC.ow <- 4.5
        
        PL$C.13CO2.ow  <- PL$C.DIC.ow * ((PL$C.DIC.ow * PL$R_ratio )/ (PL$C.DIC.ow  + (PL$C.DIC.ow * PL$R_ratio ))) # concentration HCO3 in overlying water [umol cm-3]
        PL$C.12CO2.ow  <- PL$C.DIC.ow * (PL$C.DIC.ow/(PL$C.DIC.ow  + (PL$C.DIC.ow* PL$R_ratio ))) # concentration HCO3 in overlying water [umol cm-3]
        
        
        } # DIC level
        
        if(names(S.P)[j]=="sed.rate"){PL$v<- var[k]
        PL$v.grid <- setup.prop.1D(value = PL$v, grid = PL$grid)}        # sedimentation rate
        
        
        
        #-------------------------------------------------------------------------------------------------------------------
        # Diagenetic model solution        
        #-------------------------------------------------------------------------------------------------------------------          
        output <- steady.1D(y = state, func = multi.model, parms = PL,
                            nspec = PL$N.var, positive = TRUE )
        #-------------------------------------------------------------------------------------------------------------------
        # Data storage
        #-------------------------------------------------------------------------------------------------------------------          
        # Whole sediment column d13C
        
        d13C.tot     <- (((((output$y[(7*PL$N+1):(8*PL$N)]/output$y[(6*PL$N+1):(7*PL$N)])/(PL$VPDB))-1)*1000)*PL$Frac) +                   ((PL$back.lime)*(1-PL$Frac)) 
        
        
        # Instantaneous precipitated authigenic seafloor carbonate crust, calculated from porewater DIC-C isotope composition   
        
        d13C.auth    <- ((((output$y[(1*PL$N+1):(2*PL$N)]/output$y[1:PL$N])/(PL$VPDB))-1)*1000) + PL$Delta_carb_bicarb 
        
        # Storing of diagenetic stabalized carbonate d13C data
        
        d13C.sub1   <- d13C.tot[(PL$N)]                     # at 10 meter depth
        d13C.sub2   <- (PL$Frac * median(d13C.auth[1:200]))+
          ((1-PL$Frac)* d13C.tot[(PL$N)])                     # at 0.1 meter depth
        d13C.diff1  <- c(PL$back.lime - d13C.sub1)          # at 10 meter depth, as offset (D13C) with primary carbonate d13C
        d13C.diff2  <- c(PL$back.lime - d13C.sub2)          # at 0.1 meter depth, as offset (D13C) with primary carbonate d13C
        
        if(k==1){
          Sens.sum[[t,1]] <- d13C.sub1
          Sens.sum[[t,2]] <- d13C.diff1
          Sens.sum[[t,3]] <- d13C.sub2
          Sens.sum[[t,4]] <- d13C.diff2
          Sens.sum[[t,5]] <- paste(var[k])
          Sens.sum[[t,6]] <- org.acc2[t]
          
        }else{
          Sens.sum[[(length(org.acc2)*(k-1)+t),1]] <- d13C.sub1
          Sens.sum[[(length(org.acc2)*(k-1)+t),2]] <- d13C.diff1
          Sens.sum[[(length(org.acc2)*(k-1)+t),3]] <- d13C.sub2
          Sens.sum[[(length(org.acc2)*(k-1)+t),4]] <- d13C.diff2
          Sens.sum[[(length(org.acc2)*(k-1)+t),5]] <- paste(var[k])
          Sens.sum[[(length(org.acc2)*(k-1)+t),6]] <- org.acc2[t] 
        }
        
      }
      
    }
    
    Out[[j]]<-as.data.frame(Sens.sum, stringsAsFactors = FALSE) # list collecting all sensitivity runs
    
  }
  
})


#----------------------------------------------------------------------------------------------------------------
# Results of sensitivity experiments
#----------------------------------------------------------------------------------------------------------------


Sens.F<-Out[[1]]   # sensitivity run, fraction authigenic carbonate
Sens.O2<-Out[[2]]  # sensitivity run, dissolved oxygen concentration overlying watercolumn
Sens.SO4<-Out[[3]] # sensitivity run, dissolved sulfate concentration overlying watercolumn
Sens.DIC<-Out[[4]] # sensitivity run, dissolved sulfate concentration overlying watercolumn
Sens.v<-Out[[5]]   # sensitivity run, sedimentation rate





# ============================================================================================================================
# Sensitivity tests sediment mixing
# ============================================================================================================================

PL$F.up.CH2O   <- 730.5                     # average shelf value [umol cm2- yr-1]
PL$v<- 0.2
PL$v.grid <- setup.prop.1D(value = PL$v, grid = PL$grid) 


bio.depth <- seq(0, 7.6, length.out=40)     # depth of sediment mixing


S.P<-list()
S.P$irr.0    <- c(0, 50, 200, 365.25)       # sensitivity parameter for sediment irrigation, Dale et al 2016 an van de Velde 2016
S.P$Db      <- seq(2.5, 10, length.out=4)   # sensitivity parameter for biodiffusion 



Out <- list(NA,NA,NA)                       # empty list to collect results of sensitivity tests

#-------------------------------------------------------------------------------------------------------------------
# Sensitivity test run sedimentary and oceanic parameters
#-------------------------------------------------------------------------------------------------------------------

system.time({
  
  for(j in c(1:length(S.P))){
    
    
    var<-S.P[[j]]
    
    # summary data matrix
    Sens.sum<-matrix(NA, ncol= 6, nrow= length(bio.depth)*length(var))
    colnames(Sens.sum)<-c("d13Cdia", "D13Cdia", "d13Cauth", "D13Cauth","var.sense", "Biodepth" )
    Sens.sum<-as.data.frame(Sens.sum)
    
    for(k in c(1:length(var))){
      
      for(t in c(1:length(bio.depth))){
        
        PL$x.L<- bio.depth[t] # the bioturbation depth
        
        #-------------------------------------------------------------------------------------------------------------------
        # Sensitivity of the model to the following changing variables: bio-irrigation and bio-diffusion, where the value is kept constant           between runs at zero. These latter values are considered to be representative for the post-extinction depositional setting.
        #-------------------------------------------------------------------------------------------------------------------   
        
        if(names(S.P)[j]=="irr.0"){PL$irr.0<-var[k] 
        
        PL$irr.grid <- setup.prop.1D(func = p.exp, grid = PL$grid, y.0 = PL$irr.0, y.inf = 0, x.att = PL$x.L)
        
        } else { PL$irr.0 <- 50 
        
        PL$irr.grid <- setup.prop.1D(func = p.exp, grid = PL$grid, y.0 = PL$irr.0, y.inf = 0, x.att = PL$x.L)
        
        }         # irrigation
        
        if(names(S.P)[j]=="Db"){PL$Db<-var[k]
        
        PL$Db.grid <-setup.prop.1D(func = p.exp, grid = PL$grid, y.0 = PL$Db, y.inf = 0, x.L = PL$x.L)
        
        }   # biodiffusion
        
        
        
        #-------------------------------------------------------------------------------------------------------------------
        # Diagenetic model solution        
        #-------------------------------------------------------------------------------------------------------------------          
        output <- steady.1D(y = state, func = multi.model, parms = PL,
                            nspec = PL$N.var, positive = TRUE )
        #-------------------------------------------------------------------------------------------------------------------
        # Data storage
        #-------------------------------------------------------------------------------------------------------------------          
        # Whole sediment column d13C
        
        d13C.tot     <- (((((output$y[(7*PL$N+1):(8*PL$N)]/output$y[(6*PL$N+1):(7*PL$N)])/(PL$VPDB))-1)*1000)*PL$Frac) +                   ((PL$back.lime)*(1-PL$Frac)) 
        
        
        # Instantaneous precipitated authigenic seafloor carbonate crust, calculated from porewater DIC-C isotope composition   
        
        d13C.auth    <- ((((output$y[(1*PL$N+1):(2*PL$N)]/output$y[1:PL$N])/(PL$VPDB))-1)*1000) + PL$Delta_carb_bicarb 
        
        # Storing of diagenetic stabalized carbonate d13C data
        
        d13C.sub1   <- d13C.tot[(PL$N)]                     # at 10 meter depth
        d13C.sub2   <- (PL$Frac * median(d13C.auth[1:200]))+
          ((1-PL$Frac)* d13C.tot[(PL$N)])                   # at 0.1 meter depth
        d13C.diff1  <- c(PL$back.lime - d13C.sub1)          # at 10 meter depth, as offset (D13C) with primary carbonate d13C
        d13C.diff2  <- c(PL$back.lime - d13C.sub2)          # at 0.1 meter depth, as offset (D13C) with primary carbonate d13C
        
        if(k==1){
          Sens.sum[[t,1]] <- d13C.sub1
          Sens.sum[[t,2]] <- d13C.diff1
          Sens.sum[[t,3]] <- d13C.sub2
          Sens.sum[[t,4]] <- d13C.diff2
          Sens.sum[[t,5]] <- paste(var[k])
          Sens.sum[[t,6]] <- bio.depth[t]
          
        }else{
          Sens.sum[[(length(bio.depth)*(k-1)+t),1]] <- d13C.sub1
          Sens.sum[[(length(bio.depth)*(k-1)+t),2]] <- d13C.diff1
          Sens.sum[[(length(bio.depth)*(k-1)+t),3]] <- d13C.sub2
          Sens.sum[[(length(bio.depth)*(k-1)+t),4]] <- d13C.diff2
          Sens.sum[[(length(bio.depth)*(k-1)+t),5]] <- paste(var[k])
          Sens.sum[[(length(bio.depth)*(k-1)+t),6]] <- bio.depth[t] 
        }
        
      }
      
    }
    
    Out[[j]]<-as.data.frame(Sens.sum, stringsAsFactors = FALSE) # list collecting all sensitivity runs
    
  }
  
})


#----------------------------------------------------------------------------------------------------------------
# Results of sensitivity experiments
#----------------------------------------------------------------------------------------------------------------

Sens.irr <-Out[[1]]   # sensitivity run, bio-irrigation
Sens.Db  <-Out[[2]]   # sensitivity run, bio-diffusion




#----------------------------------------------------------------------------------------------------------------
# Plots of combined sensitivity experiments
#----------------------------------------------------------------------------------------------------------------

var<-c("D13Cdia", "D13Cauth") # end-member carbonate d13C values of interest

for(i in c(1:2)){
  if(var[i] == "D13Cdia"){
    lims <- c(-3,3)}else{ 
      lims <- c(-5,10)}
  
  # Plot of authigenic carbonate fraction sensitivity experiment
  
  A<-ggplot(Sens.F, aes(x=FOC, y=Sens.F[,var[i]], group=var.sense, color=as.factor(var.sense) ))+geom_line()+
    scale_color_discrete(name="fraction authigenic \n carbonate")+
    ylab(expression(atop(paste(Delta^13*C[primary-bulk]~"(VPDB)"))))+
    xlab(expression(atop(paste("F"[OC]~(mu*mol~cm^-3~y^-1)))))+
    ggtitle("(a)")+
    ylim(lims)+
    theme(
      axis.line.x = element_line(color = "black"), 
      axis.line.y = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = c(0.75,0.75),
      legend.title = element_text(size=7),
      legend.text = element_text(size=7),
      plot.title = element_text(hjust=(0)))+
    guides(colour=guide_legend(ncol=2))
  
  # Plot of oxygen level sensitivity experiment
  
  
  B<-ggplot(Sens.O2, aes(x=FOC, y=Sens.O2[,var[i]], group=var.sense, color=as.factor(var.sense) ))+geom_line()+
    scale_color_discrete(name=expression(paste(O[2]~(mu*mol~cm^-3))))+
    ylab(expression(atop(paste(Delta^13*C[primary-bulk]~"(VPDB)"))))+
    xlab(expression(atop(paste("F"[OC]~(mu*mol~cm^-3~y^-1)))))+
    ggtitle("(b)")+
    ylim(lims)+
    theme(
      axis.line.x = element_line(color = "black"), 
      axis.line.y = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = c(0.75,0.75),
      legend.title = element_text(size=7),
      legend.text = element_text(size=7),
      plot.title = element_text(hjust=(0)))+
    guides(colour=guide_legend(ncol=2))
  
  # Plot of sulfate level sensitivity experiment
  
  C<-ggplot(Sens.SO4, aes(x=FOC, y=Sens.SO4[,var[i]], group=var.sense, color=as.factor(var.sense) ))+geom_line()+
    scale_color_discrete(name=expression(paste(SO[4]~(mu*mol~cm^-3))))+
    ylab(expression(atop(paste(Delta^13*C[primary-bulk]~"(VPDB)"))))+
    xlab(expression(atop(paste("F"[OC]~(mu*mol~cm^-3~y^-1)))))+
    ggtitle("(c)")+
    ylim(lims)+
    theme(
      axis.line.x = element_line(color = "black"), 
      axis.line.y = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = c(0.75,0.75),
      legend.title = element_text(size=7),
      legend.text = element_text(size=7),
      plot.title = element_text(hjust=(0)))+
    guides(colour=guide_legend(ncol=2))
  
  # Plot of DIC level sensitivity experiment
  
  D<-ggplot(Sens.DIC, aes(x=FOC, y=Sens.DIC[,var[i]], group=var.sense, color=as.factor(var.sense) ))+geom_line()+
    scale_color_discrete(name=expression(paste(DIC~(mu*mol~cm^-3))))+
    ylab(expression(atop(paste(Delta^13*C[primary-bulk]~"(VPDB)"))))+
    xlab(expression(atop(paste("F"[OC]~(mu*mol~cm^-3~y^-1)))))+
    ggtitle("(d)")+
    ylim(lims)+
    theme(
      axis.line.x = element_line(color = "black"), 
      axis.line.y = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = c(0.75,0.75),
      legend.title = element_text(size=7),
      legend.text = element_text(size=7),
      plot.title = element_text(hjust=(0)))+
    guides(colour=guide_legend(ncol=2))
  
  
  
  
  # Plot of sedimentation rate sensitivity experiment
  
  E<-ggplot(Sens.v, aes(x=FOC, y=Sens.v[,var[i]], group=var.sense, color=as.factor(var.sense) ))+geom_line()+
    scale_color_discrete(name=expression(paste("sedimentation rate"~(cm~y^-1))))+
    ylab(expression(atop(paste(Delta^13*C[primary-bulk]~"(VPDB)"))))+
    xlab(expression(atop(paste("F"[OC]~(mu*mol~cm^-3~y^-1)))))+
    ggtitle("(e)")+
    ylim(lims)+
    theme(
      axis.line.x = element_line(color = "black"), 
      axis.line.y = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = c(0.7,0.75),
      legend.title = element_text(size=7),
      legend.text = element_text(size=7),
      plot.title = element_text(hjust=(0)))+
    guides(colour=guide_legend(ncol=2))
  
  F<-ggplot(Sens.irr, aes(x=Biodepth, y=Sens.irr[,var[i]], group=var.sense, color=as.factor(var.sense) ))+geom_line()+
    scale_color_discrete(name=expression(paste("Bio-irrigation"~(y^-1))))+
    ylab(expression(atop(paste(Delta^13*C[primary-bulk]~"(VPDB)"))))+
    xlab("attenuation depth coefficient (cm)")+
    
    # Modern Van de Velde and Meysman 2016 Aqua Chem
    geom_vline(xintercept=2)+ 
    annotate("text", label="Modern", y=0,  x=2.1 , size=3.5, angle=90)+
    
    # Paleozoic Dale et al 2016 GCA
    geom_vline(xintercept=1, linetype = 2)+ 
    annotate("text", label="Pre-extinction", y=0,  x=1.1 , size=3.5, angle=90)+
    
    # Earliest Paleozoic (no burrowers) Dale et al 2016 GCA 
    geom_vline(xintercept=0, linetype = 3)+ 
    annotate("text", label="Post_extinction", y=0,  x=0.1 , size=3.5, angle=90)+
    
    ggtitle("(f)")+
    
    ylim(lims)+
    theme(
      axis.line.x = element_line(color = "black"), 
      axis.line.y = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = c(0.75,0.75),
      legend.title = element_text(size=7),
      legend.text = element_text(size=7),
      plot.title = element_text(hjust=(0)))+
    guides(colour=guide_legend(ncol=2))
  
  
  # Plot of biodiffusion
  
  
  G<-ggplot(Sens.Db, aes(x=Biodepth, y=Sens.Db[,var[i]], group=var.sense, color=as.factor(var.sense) ))+geom_line()+
    scale_color_discrete(name=expression(paste("Bio-diffusion"~(cm^2~y^-1))))+
    ylab(expression(atop(paste(Delta^13*C[primary-bulk]~"(VPDB)"))))+
    xlab("mixing depth (cm)")+
    
    # earliest Palaeozoic Dale et al 2016 GCA
    geom_vline(xintercept=0, linetype = 4)+
    annotate("text", label="Post-extinction", y=0,  x=0.1 , size=3.5, angle=90)+
    
    # Palaeozoic Dale et al 2016 GCA
    geom_vline(xintercept=2, linetype = 3)+
    annotate("text", label="Pre-extinction", y=0,  x=2.1 , size=3.5, angle=90)+
    
    # Weak mixing vandevelde & Meysman 2016 Aqua Chem
    geom_vline(xintercept=3.6, linetype = 2)+ 
    annotate("text", label="weak mixing", y=0,  x=3.7 , size=3.5, angle=90)+
    
    # Strong mixing vandevelde & Meysman 2016
    geom_vline(xintercept=7.6, linetype = 1)+ 
    annotate("text", label="strong mixing", y=0,  x=7.7 , size=3.5, angle=90)+
    
    ggtitle("(g)")+
    ylim(lims)+
    theme(
      axis.line.x = element_line(color = "black"), 
      axis.line.y = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = c(0.75,0.75),
      legend.title = element_text(size=7),
      legend.text = element_text(size=7),
      plot.title = element_text(hjust=(0)))+
    guides(colour=guide_legend(ncol=2))
  
  grid.arrange(A,B,C,D,E,F,G)
  
  FOC_sens<-arrangeGrob(A,B,C,D,E,F,G)
  
  if(var[i]=="D13Cdia"){
    ggsave("Figure5DIC.pdf",FOC_sens, height= 24, width=24, units="cm" )} # saves produced graphs of solely ongoing diagenesis with depth
  
  if(var[i]=="D13Cauth"){
    ggsave("Figure6DIC.pdf",FOC_sens, height= 24, width=24, units="cm" )} # saves produced graphs of solely seafloor carbonate crust 
  
}






#=================================================================================================================
# Timeseries analysis
#=================================================================================================================


#------------------------------------------------------------------------------------------------------------------------- 
# Reset parameters to Permian conditions
#------------------------------------------------------------------------------------------------------------------------- 

#DIC
PL$C.DIC.ow <- 4.5
PL$C.13CO2.ow  <- PL$C.DIC.ow * ((PL$C.DIC.ow * PL$R_ratio )/ (PL$C.DIC.ow  + (PL$C.DIC.ow * PL$R_ratio ))) # concentration HCO3 in overlying water [umol cm-3]
PL$C.12CO2.ow  <- PL$C.DIC.ow * (PL$C.DIC.ow/(PL$C.DIC.ow  + (PL$C.DIC.ow* PL$R_ratio ))) # concentration HCO3 in overlying water [umol cm-3]

# SO4
PL$C.SO4.ow <- 4

# O2
PL$C.O2.ow <- 0.28

# Diagenetic/Authigenic fraction
PL$Frac <- 0.2

# Bio-diffusion
PL$x.L <- 2 # Phanerozoic mixing depth [cm] Dale et al. 2016 GCA
PL$Db <- 5
PL$Db.grid <-setup.prop.1D(func = p.exp, grid = PL$grid, y.0 = PL$Db, y.inf = 0, x.L = PL$x.L)

# Bio-irrigation
PL$z.bio <- 1 # Phanerozoic attenuation depth coefficient [cm] Dale et al. 2016 GCA
PL$irr.0 <- 50
PL$irr.grid <- setup.prop.1D(func = p.exp, grid = PL$grid, y.0 = PL$irr.0, y.inf = 0, x.att = PL$z.bio)

# Sedimentation rate
PL$v<- 0.2
PL$v.grid <- setup.prop.1D(value = PL$v, grid = PL$grid)

#------------------------------------------------------------------------------------------------------------------------- 

NPop<-50 # sample population size for generating spatial heterogeneous OC accumulation

# Parameters for timeseries taken from the actual d13C data of the P-Tr rock sequences in Iran and China (CarbTrends.R)

#-------------------------------------------------------------------------------
# Iran
bulk.lime <- ModelInput.Iran$M      # bulk-rock carebon isotope value, Iran
IQR.seq   <- ModelInput.Iran$IQR.M  # carbon isotope IQR, Iran

#-------------------------------------------------------------------------------
# China
#bulk.lime <- ModelInput.China$M      # bulk-rock carebon isotope value, China
#IQR.seq   <- ModelInput.China$IQR.M  # carbon isotope IQR, China

#-------------------------------------------------------------------------------

# Linear relation between dispersion and modulation of OC flux size

up.org.f  <- 4 # high global OC flux Induan, after Algeo et al., 2013, given as factor
low.org.f <- 1 # low global OC flux Permian, after Algeo et al., 2013, given as factor

factor  <- ((up.org.f-low.org.f)/(max(IQR.seq)-min(IQR.seq)))*IQR.seq +(low.org.f-(((up.org.f-low.org.f)/(max(IQR.seq )-min(IQR.seq)))*min(IQR.seq)))

# Inverse linear relation between spatial dispersion and modulation of OC flux density distribution (i.e. spatial-heterogeneity of OC accumulation)

up.diss.f  <- 5 # high OC dispersion factor
low.diss.f <- 2 # low OC dispersion factor

factor2 <- ((up.diss.f-low.diss.f)/(max(factor)-min(factor)))*factor+(low.diss.f-(((up.diss.f-low.diss.f)/(max(factor)-min(factor)))*min(factor)))

# Linear relation between spatial dispersion of OC and authigenic seafloor carbonate crust formation

up.rho  <- 0.99 # upper probability of not forming authigenic seafloor carbonate crust formation
low.rho <- 0.90  # lower probability of not forming authigenic seafloor carbonate crust formation

factor3 <- ((low.rho-up.rho)/(5-3))*factor2+(up.rho-(((low.rho-up.rho)/(5-3))*(3)))
factor3[factor3>1]<- 1

# OC flux calculated based on the sensitivity of carbon isotope composition towards changing OC fluxes and a 4 mmol kg- SO4 level 

org.acc<-(((Sens.v[which(Sens.v$var.sense==0.2),]))[which((Sens.v[which(Sens.v$var.sense==0.2),])$D13Cdia == max(((Sens.v[which(Sens.v$var.sense==0.2),])$D13Cdia))
),] )$FOC 

org.acc<-(((Sens.SO4[which(Sens.SO4$var.sense==4),]))                      
          [which((Sens.SO4[which(Sens.SO4$var.sense==4),])$D13Cdia == max(((Sens.SO4[which(Sens.SO4$var.sense==4),])$D13Cdia))
          ),] )$FOC 

# Empty vectors

d13C.Dia  <-rep(NA, length.out=c(NPop*length(factor2)))
d13C.Auth <-rep(NA, length.out=c(NPop*length(factor2)))
d13C.Bulk <-rep(NA, length.out=c(NPop*length(factor2)))
d13C.sw   <-rep(NA, length.out=c(NPop*length(factor2)))
FOC       <-rep(NA, length.out=c(NPop*length(factor2)))


#=================================================================================================================
# Timeseries model run
#=================================================================================================================

system.time({
  
  
  for(t in c(1:length(factor2))){
    
    # print progress
    if(t %% 10==0) { cat(paste0("iteration: ", i, "\n")) }  
    
    
    d13C.sub1 <-rep(NA, length.out= NPop)# vectors for storing generated sampled sediment carbonate del13C
    d13C.sub2 <-rep(NA, length.out= NPop)
    d13C.sub3 <-rep(NA, length.out= NPop)
    d13C.sub4 <-rep(NA, length.out= NPop)
    d13C.sub5 <-rep(NA, length.out= NPop)
    
    d13C.tot<-rep(NA, length.out= NPop*PL$N)
    d13C.auth<-rep(NA, length.out= NPop*PL$N)
    
    for(i in c(1:NPop)){
      
      #-------------------------------------------------------------------------------------------------------------------------     
      # The OC flux and spatial dispersion
      #-------------------------------------------------------------------------------------------------------------------------        
      
      mu <- log(org.acc*factor[t])       # log-transformed mean OC flux
      sigma <- log(factor2[t])           # log-transformed standard devation of OC flux
      
      set.seed(i*t)
      PL$F.up.CH2O <- rlnorm(1,mu,sigma) # random generated OC flux from lognormal density distribution
      
      #------------------------------------------------------------------------------------------------------------------------- 
      
      #-------------------------------------------------------------------------------------------------------------------------  
      # Extinction of benthic fauna at the exinction horizon and its associated effect on bio-irrigation and bio-diffusion 
      #-------------------------------------------------------------------------------------------------------------------------       
      if(t > 60){
        PL$irr.0 <- 0               # sediment bio-irrigation
        PL$irr.grid <- setup.prop.1D(func = p.exp, grid = PL$grid, y.0 = PL$irr.0, y.inf = 0, x.att = PL$x.L)
        PL$Db <- 0                  # sediment bio-difussion
        PL$Db.grid <- setup.prop.1D(func = p.exp, grid = PL$grid, y.0 = PL$Db, y.inf = 0, x.L = PL$x.L)}
      #-------------------------------------------------------------------------------------------------------------------------  
      
      #-------------------------------------------------------------------------------------------------------------------------  
      # Sedimentation rate change over the Permian-Triassic boundary
      #-------------------------------------------------------------------------------------------------------------------------  
      
      if(t > 70){
        PL$v<- 0.4
        PL$v.grid <- setup.prop.1D(value = PL$v, grid = PL$grid)}
      
      
      #------------------------------------------------------------------------------------------------------------------------- 
      # The sold and fluid carbonate carbon isotope chemistry
      #-------------------------------------------------------------------------------------------------------------------------         
      
      # Correcting the observed bulk rock end-member isotope signal with the isotope offset obtained by the sensitivity 
      # Experiments with changing sedimentation rates
      
      PL$back.lime <- bulk.lime[t] + 
        (Sens.v[which(Sens.v$var.sense == PL$v & 
                        (abs(Sens.v$FOC-(org.acc*factor[t]))==
                           min(abs(Sens.v$FOC-(org.acc*factor[t]))))),])$D13Cdia # primary solid carbonate d13C 
      
      PL$back.fluid   <- PL$back.lime-PL$Delta_carb_bicarb # primary fluid bicarbonate d13C 
      
      # 13C abundances of dissolved HCO and fixed in carbonate mineral lattice by 20 degree Celcius
      
      PL$Rcarb_ratio <- PL$VPDB * (PL$back.lime/1000+1)    # precipitated carbonate carbon isotope ratio Permian ocean water
      PL$R_ratio     <- PL$VPDB * (1+(PL$back.fluid/1000)) # bicarbonate carbon isotope ratio Permian ocean water
      
      # Calciumcarbonate d13C correction
      
      PL$F.up.12Carb <-PL$F.up.Carb * (1/(1 + PL$Rcarb_ratio))               # concentration 12CO3 in  carbonate [umol cm2- yr-1]
      PL$F.up.13Carb <-PL$F.up.Carb * (PL$Rcarb_ratio/(1 + PL$Rcarb_ratio))  # concentration 13CO3 in  carbonate [umol cm2- yr-1]
      
      # Water d13C correction
      
      PL$C.12CO2.ow  <- PL$C.DIC.ow * (PL$C.DIC.ow/
                                      (PL$C.DIC.ow  + (PL$C.DIC.ow* PL$R_ratio )))# concentration 12HCO3 in overlying water [umol cm-3]
      PL$C.13CO2.ow  <- PL$C.DIC.ow * ((PL$C.DIC.ow * PL$R_ratio )/ 
                                      (PL$C.DIC.ow  + (PL$C.DIC.ow * PL$R_ratio )))# concentration 13HCO3 in overlying water [umol cm-3]
      
      
      #-------------------------------------------------------------------------------------------------------------------------------
      output <- steady.1D(y = state, func = multi.model, parms = PL,
                          nspec = PL$N.var, positive = TRUE )
      #-------------------------------------------------------------------------------------------------------------------------------
      
      #-------------------------------------------------------------------------------------------------------------------------------   
      # Data storage   
      #-------------------------------------------------------------------------------------------------------------------------------    
      if(i==1){# the first cycle of the timeseries
        
        # Whole sediment column d13C 
        
        d13C.tot[1:(PL$N*i)] <-(((((output$y[(7*PL$N+1):(8*PL$N)]/output$y[(6*PL$N+1):(7*PL$N)])/
                                     (PL$VPDB))-1)*1000)*PL$Frac) + ((PL$back.lime)*(1-PL$Frac)) 
        
        # Instantaneous precipitated authigenic seafloor carbonate crust, calculated from DIC C isotope composition
        
        d13C.auth[1:(PL$N*i)] <- ((((output$y[(1*PL$N+1):(2*PL$N)]/
                                       output$y[1:PL$N])/(PL$VPDB))-1)*1000) + PL$Delta_carb_bicarb 
        
        d13C.sub1[i]   <- d13C.tot[(PL$N)]        # at 10 meter depth
        d13C.sub2[i]   <- median(d13C.auth[1:200]) # at 0.1 meter depth
        
        # Assigning a probability to authigenic carbonate formation and its impact on bulk rock carbon isotope composition   
        
        if( sample(x=c(0, 1), size=1 , prob=c(factor3[t],1-factor3[t]))==1){
          d13C.sub3[i]   <- 0.8 * d13C.sub1[i] + 0.2*d13C.sub2[i]}else{
            d13C.sub3[i]   <- d13C.sub1[i]  
          }  
        
        d13C.sub4[i]   <- PL$back.fluid # the d13C seawater value
        d13C.sub5[i]   <- PL$F.up.CH2O  # the average OC flux value
        
      }else{# the consecutive cycles of the timeseries
        
        # Whole sediment column d13C     
        
        d13C.tot[(((PL$N*i)-PL$N)+1):(PL$N*i)]<-(((((output$y[(7*PL$N+1):(8*PL$N)]/
                                                       output$y[(6*PL$N+1):(7*PL$N)])/(PL$VPDB))-1)*1000)*PL$Frac)+ 
          ((PL$back.lime)*(1-PL$Frac)) 
        
        # Instantaneous precipitated authigenic seafloor carbonate crust, calculated from DIC C isotope composition    
        
        d13C.auth[(((PL$N*i)-PL$N)+1):(PL$N*i)]<- ((((output$y[(1*PL$N+1):(2*PL$N)]/
                                                        output$y[1:PL$N])/(PL$VPDB))-1)*1000) + PL$Delta_carb_bicarb 
        
        d13C.sub1[i]   <- d13C.tot[(PL$N)*i]                                                # at 10 meter depth
        d13C.sub2[i]   <- median(d13C.auth[(((PL$N)*i)-((PL$N))):(((PL$N)*i)-((PL$N)-20))]) # at 0.1 meter depth
        
        # Assigning a probability to authigenic carbonate formation and its impact on bulk rock carbon isotope composition         
        
        if( sample(x=c(0, 1), size=1 , prob=c(factor3[t],1-factor3[t]))==1){
          d13C.sub3[i]   <- 0.8 * d13C.sub1[i] + 0.2*d13C.sub2[i]}else{
            d13C.sub3[i]   <- d13C.sub1[i]  
          } 
        
        d13C.sub4[i]   <- PL$back.fluid # the d13C seawater value
        d13C.sub5[i]   <- PL$F.up.CH2O  # the average OC flux value
      }
    }
    
    if(t==1){
      d13C.Dia[t:NPop]  <- d13C.sub1
      d13C.Auth[t:NPop] <- d13C.sub2
      d13C.Bulk[t:NPop] <- d13C.sub3
      d13C.sw[t:NPop]   <- d13C.sub4
      FOC[t:NPop]       <- d13C.sub5 }else{
        
        d13C.Dia[c(NPop*t+1-NPop):c(NPop*t)]  <- d13C.sub1
        d13C.Auth[c(NPop*t+1-NPop):c(NPop*t)] <- d13C.sub2
        d13C.Bulk[c(NPop*t+1-NPop):c(NPop*t)] <- d13C.sub3
        d13C.sw[c(NPop*t+1-NPop):c(NPop*t)]   <- d13C.sub4
        FOC[c(NPop*t+1-NPop):c(NPop*t)]       <- d13C.sub5 
        
      }
    
  }})

#=================================================================================================================
#=================================================================================================================

#----------------------------------------------------------------------------------------------------------------------------------
# Data conversion for plotting
#----------------------------------------------------------------------------------------------------------------------------------
TimeSeries.com.Iran <- data.frame(d13C.Bulk, grid = rep(ModelInput.Iran$x,each=NPop))

TimeSeries.sw.Iran <- data.frame(d13C.sw=d13C.sw[seq(1,10100, NPop)], grid = ModelInput.Iran$x) # seawater DIC d13C

IPR<-function(x, seR=seq(0,1,0.025)){(quantile(x,seR, na.rm=TRUE))[40]-(quantile(x,seR, na.rm=TRUE))[2]} # interpercentile range [95%]
IQR<-function(x, seR=seq(0,1,0.25)){(quantile(x,seR, na.rm=TRUE))[4]-(quantile(x,seR, na.rm=TRUE))[2]}   # interquartile range [50%]

# Median trend lines, functions obtained from CarbTrends.R script
BootMedian.Iran<-BootSlide(n=999,data=TimeSeries.com.Iran$d13C.Bulk ,time=TimeSeries.com.Iran$grid, grid=grid.Iran, windowsize=windowsize, sampsize=10, method=median)
BootIQR.Iran<-BootSlide(n=999,data=TimeSeries.com.Iran$d13C.Bulk ,time=TimeSeries.com.Iran$grid, grid=grid.Iran, windowsize=windowsize, sampsize=10, method=IQR)
BootIQR_95.Iran<-BootSlide(n=999,data=TimeSeries.com.Iran$d13C.Bulk ,time=TimeSeries.com.Iran$grid, grid=grid.Iran, windowsize=windowsize, sampsize=10, method=IPR)

BootSum.Iran<-cbind(rbind(BootMedian.Iran, BootIQR.Iran, BootIQR_95.Iran), c(rep("Median", nrow(BootMedian.Iran)), rep("IQR (50 %)", nrow(BootIQR.Iran)), rep("IQR (95%)", nrow(BootIQR_95.Iran) ))) # summary

# Weighted color, functions obtained from CarbTrends.R script
colorMed.Iran<-colorWeight(data=BootMedian.Iran, grid=grid.Iran)
colorIQR.Iran<-colorWeight(data=BootIQR.Iran, grid=grid.Iran)
colorIQR_95.Iran<-colorWeight(data=BootIQR_95.Iran, grid=grid.Iran)


CI.sum.Iran<-cbind(rbind((colorMed.Iran[[1]]), (colorIQR.Iran[[1]]), (colorIQR_95.Iran[[1]])), label=c(rep("Median", nrow((colorMed.Iran[[1]]))), rep("IQR (50 %)", nrow((colorIQR.Iran[[1]]))), rep("IPR (95%)", nrow((colorIQR_95.Iran[[1]])) )), stringsAsFactors=FALSE)

color.sum.Iran<-cbind(rbind((colorMed.Iran[[2]]), (colorIQR.Iran[[2]]), (colorIQR_95.Iran[[2]])), label=c(rep("Median", nrow((colorMed.Iran[[2]]))), rep("IQR (50 %)", nrow((colorIQR.Iran[[2]]))), rep("IPR (95%)", nrow((colorIQR_95.Iran[[2]])))), stringsAsFactors=FALSE)

#----------------------------------------------------------------------------------------------------------------------------------

palette<-colorRampPalette(c("#a6bddb", "#034e7b"), bias=2)(20) # color palette

# Plot of model results

ggplot(CI.sum.Iran, aes(x=x, y=M,  alpha=w3^3))+
  
  facet_grid(label~., scales="free_x")+
  
  # Annotation for the extinction horizon
  geom_vline(xintercept=(6), lty=2)+
  
  # Annotation for the Permian-Triassic boundary
  geom_vline(xintercept=7, lty=1)+
  
  # Smoothed trendline
  geom_tile(data=color.sum.Iran, aes(x=x, y=y, fill=dens.scaled, alpha=alpha.factor))+
  scale_fill_gradientn("dens.scaled", colours=palette, guide=FALSE)+ 
  scale_alpha_continuous(range=c(0.001, 1), guide=FALSE)+
  geom_line(color="white")+
  
  # Axis  
  
  ylab(expression(atop(paste(delta^13*C~"(\211 VPDB)"))))+
  scale_y_continuous(limits=c(-10, 10), expand = c(0, 0))+
  
  # Overwrites the labels for the aes x to create the names of the biozones 
  
  scale_x_discrete("biozones", limits=seq(0, 11,1) , expand = c(0, 0), labels=c("",paste(LETTERS[1:11])))+ 
  coord_cartesian(xlim=seq(0, 11,1))+
  
  # Chronological boxes for the systems
  
  geom_rect(aes(xmax=bounds.system[1], xmin=0, ymax=-9, ymin=-10), fill=strat.colors.system[2])+
  geom_rect(aes(xmax=bounds.system[2], xmin=bounds.system[1], ymax=-9, ymin=-10), fill=strat.colors.system[1])+
  
  
  
  # Chronological boxes for the stages
  
  geom_rect(aes(xmax=bounds.stage[1], xmin=0,  ymax=-8, ymin=-9), fill=strat.colors.stage[3])+
  geom_rect(aes(xmax=bounds.stage[2], xmin=bounds.stage[1],  ymax=-8, ymin=-9), fill=strat.colors.stage[2])+
  geom_rect(aes(xmax=bounds.stage[3], xmin=bounds.stage[2],  ymax=-8, ymin=-9), fill=strat.colors.stage[1])+
  
  
  theme(
    axis.line.x = element_line(color = "black"), 
    axis.line.y = element_line(color = "black"),
    legend.title = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    legend.text = element_text(size=5), 
    legend.position = c(0.18,0.29), 
    legend.key.size = unit(3,"mm"),
    legend.background = element_blank(),
    
    axis.text.x  = element_text(face = "bold", hjust = 2.5, size = 5), 
    axis.title.x  = element_text(size = 7), 
    axis.text.y = element_text(size=7), 
    axis.title.y = element_text(size =7),
    plot.title = element_text(lineheight=.8, face="bold")
    
  )

ggsave("Figure8.pdf", height= 6, width=4, units="in") # saving model graph

#---------------------------------------------------------------------------
#China, Figure Extra

TimeSeries.com.China <- data.frame(d13C.Bulk, grid = rep(ModelInput.China$x,each=NPop))

TimeSeries.sw.China <- data.frame(d13C.sw=d13C.sw[seq(1,8100, 100)], grid = ModelInput.China$x) # seawater DIC d13C



# Median trend lines, functions obtained from CarbTrends.R script
BootMedian.China<-BootSlide(n=999,data=TimeSeries.com.China$d13C.Bulk ,time=TimeSeries.com.China$grid, grid=ModelInput.China$x, windowsize=windowsize, sampsize=10, method=median)
BootIQR.China<-BootSlide(n=999,data=TimeSeries.com.China$d13C.Bulk ,time=TimeSeries.com.China$grid, grid=ModelInput.China$x, windowsize=windowsize, sampsize=10, method=IQR)
BootIQR_95.China<-BootSlide(n=999,data=TimeSeries.com.China$d13C.Bulk ,time=TimeSeries.com.China$grid, grid=ModelInput.China$x, windowsize=windowsize, sampsize=10, method=IPR)

BootSum.China<-cbind(rbind(BootMedian.China, BootIQR.China, BootIQR_95.China), c(rep("Median", nrow(BootMedian.China)), rep("IQR (50 %)", nrow(BootIQR.China)), rep("IQR (95%)", nrow(BootIQR_95.China) ))) # summary

# Weighted color, functions obtained from CarbTrends.R script

colorMed.China<-colorWeight(data=BootMedian.China, grid=ModelInput.China$x)
colorIQR.China<-colorWeight(data=BootIQR.China, grid=ModelInput.China$x)
colorIQR_95.China<-colorWeight(data=BootIQR_95.China, grid=ModelInput.China$x)


CI.sum.China<-cbind(rbind((colorMed.China[[1]]), (colorIQR.China[[1]]), (colorIQR_95.China[[1]])), label=c(rep("Median", nrow((colorMed.China[[1]]))), rep("IQR (50 %)", nrow((colorIQR.China[[1]]))), rep("IPR (95%)", nrow((colorIQR_95.China[[1]])) )), stringsAsFactors=FALSE)

color.sum.China<-cbind(rbind((colorMed.China[[2]]), (colorIQR.China[[2]]), (colorIQR_95.China[[2]])), label=c(rep("Median", nrow((colorMed.China[[2]]))), rep("IQR (50 %)", nrow((colorIQR.China[[2]]))), rep("IPR (95%)", nrow((colorIQR_95.China[[2]])))), stringsAsFactors=FALSE)

#----------------------------------------------------------------------------------------------------------------------------------

palette<-colorRampPalette(c("#a6bddb", "#034e7b"), bias=2)(20) # color palette

# plot of model results

ggplot(CI.sum.China, aes(x=x, y=M,  alpha=w3^3))+
  
  facet_grid(label~., scales="free_x")+
  
  # annotation for the extinction horizon
  geom_vline(xintercept=(if(subject == "China"){6} else {5.9}), lty=2)+
  
  # annotation for the Permian-Triassic boundary
  geom_vline(xintercept=7, lty=1)+
  
  # smoothed trendline
  geom_tile(data=color.sum.China, aes(x=x, y=y, fill=dens.scaled, alpha=alpha.factor))+
  scale_fill_gradientn("dens.scaled", colours=palette, guide=FALSE)+ 
  scale_alpha_continuous(range=c(0.001, 1), guide=FALSE)+
  geom_line(color="white")+
  
  # Axis  
  
  ylab(expression(atop(paste(delta^13*C~"(\211 VPDB)"))))+
  scale_y_continuous(limits=c(-10, 10), expand = c(0, 0))+
  
  # Overwrites the labels for the aes x to create the names of the biozones 
  
  scale_x_discrete("biozones", limits=seq(0, 11,1) , expand = c(0, 0), labels=c("",paste(LETTERS[1:11])))+ 
  coord_cartesian(xlim=seq(0, 11,1))+
  
  # Chronological boxes for the systems
  
  geom_rect(aes(xmax=bounds.system[1], xmin=0, ymax=-9, ymin=-10), fill=strat.colors.system[2])+
  geom_rect(aes(xmax=bounds.system[2], xmin=bounds.system[1], ymax=-9, ymin=-10), fill=strat.colors.system[1])+
  
  
  
  # Chronological boxes for the stages
  
  geom_rect(aes(xmax=bounds.stage[1], xmin=0,  ymax=-8, ymin=-9), fill=strat.colors.stage[3])+
  geom_rect(aes(xmax=bounds.stage[2], xmin=bounds.stage[1],  ymax=-8, ymin=-9), fill=strat.colors.stage[2])+
  geom_rect(aes(xmax=bounds.stage[3], xmin=bounds.stage[2],  ymax=-8, ymin=-9), fill=strat.colors.stage[1])+
  
  
  theme(
    axis.line.x = element_line(color = "black"), 
    axis.line.y = element_line(color = "black"),
    legend.title = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x=element_blank(),
    legend.text = element_text(size=5), 
    legend.position = c(0.18,0.29), 
    legend.key.size = unit(3,"mm"),
    legend.background = element_blank(),
    
    axis.text.x  = element_text(face = "bold", hjust = 2.5, size = 5), 
    axis.title.x  = element_text(size = 7), 
    axis.text.y = element_text(size=7), 
    axis.title.y = element_text(size =7),
    plot.title = element_text(lineheight=.8, face="bold")
    
  )

ggsave("FigureExtra.pdf", height= 6, width=4, units="in") # saving model graph

################################################################################################################################
################################################################################################################################
################################################################################################################################
