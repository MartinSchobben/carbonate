##################################################################################################################
# Data visualization and statistics on the Iranian and Chinese del13C Permian-Triassic stratigraphic profiles 
# Construction of Figures 1, 3 and 4 is followed by Figure 2 and are presented in the work titled; Latest Permian stable carbon isotope variability traces heterogeneous organic carbon accumulation and authigenic carbonate formation, published in Climate of the Past Discussions. This R script is complemented by a diagenetic model named; CarDiaModel.R, which is required to construct the Figures 5, 6 and 8 of the same work.
######################################################################################################################

# required packages
library(ggplot2)    # for data plotting
library(plyr)       # for split array-apply functions
library(reshape2)   # for melting of a dataframe
library(gridExtra)  # for multi-graph option
library(grid) # grid.draw-function


# =============================================================================
# Global settings
# =============================================================================


PTr<-read.csv("Schobben_SupplementaryData2.csv")
# selecting d13C values with biostratigraphic control 
Good<-complete.cases(PTr$del13C,PTr$biozone.label) 
PTr<-PTr[Good,]


# =============================================================================
# Construction of a dimensionless timeframe
# =============================================================================


# biozones
A<-c(-4322, NA, -3000, NA, NA, NA, -4000, NA, 1.00)
B<-c(-1622, NA, -1870, NA, -2170, NA, -1403, NA, 0.50)
C<-c(-782, NA, -1442, NA, -1598, NA, -579, -4101, 0.35)  
D<-c(-482, NA, -819, -470, -973, -1067, NA, -2581, 0.86)
E<-c(-388, -398, -389, -190, -303, -287, -297, -303, 0.86)
F<-c(-60, -22, -31, NA, -46, -42, -47, -30, 0.15)
G<-c(0,    0,   0,  0,   0,   0,   0,    4, 0.05)
H<-c(138,  185, NA, 145, 132, 210, 159,  28, 0.02)
I<-c(237,  371, NA, 1025, 546, 715, 443,  40, 0.1)
J<-c(878,  NA, NA, NA, NA, NA, NA,  913, 0.58)
K<-c(1438,  NA, NA, NA, NA, NA, NA,  4113, 0.30)  

CONO<-rbind(A, B, C, D, E, F, G, H, I, J, K) # biozone lower boundaries, cm relative to the base of the Boundary Clay, encountered in P-Tr profiles situated in China as well as Iran 

colnames(CONO)<-c(unique(as.character(PTr$locality)), "Duration") # data frame consisting of the FAD and LAD of conodont specimens marking conodont zones, here represented as capital letters. The boundary represents the lower (oldest) boundary of the biozone. Last column refers

DUR<-CONO[,9] # biozone durations [My]
CONO<-CONO[,-9] # biozone stratigraphic height [cm]
colnames(CONO)<-unique(PTr$locality)

THICK<-rbind(CONO[c(2:nrow(CONO)),]-CONO[c(1:(nrow(CONO)-1)),], rep(NA, ncol(CONO))) # biozone stratigraphic thickness [cm]
colnames(THICK)<-unique(PTr$locality)
rownames(THICK)<-c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K")


dim_height<-c()

for(i in c(1:nrow(PTr))){# 
  
  dim_height<-c(dim_height,
                ((PTr$height[i]-CONO[as.character(PTr$biozone.label[i]),as.character(PTr$locality[i])]) / (THICK[as.character(PTr$biozone.label[i]),as.character(PTr$locality[i])])))
  
}# conversion of stratigraphic heights to dimensionless heights. The dimensionless height is assigned to the d13C values based on the relative distance towards the lower boundary of the corresponding biostratigraphic unit.

PTr<-cbind(PTr, dim_height)

#if(PTr$code[i] %in% rownames(CONO))

multiples <- # dimensionless heights additive
  outer(PTr$biozone.label=="A", 0) +
  outer(PTr$biozone.label=="B", 1) +
  outer(PTr$biozone.label=="C", 2) +
  outer(PTr$biozone.label=="D", 3) +
  outer(PTr$biozone.label=="E", 4) +
  outer(PTr$biozone.label=="F", 5) +
  outer(PTr$biozone.label=="G", 6) +
  outer(PTr$biozone.label=="H", 7) +
  outer(PTr$biozone.label=="I", 8) +
  outer(PTr$biozone.label=="J", 9) +
  outer(PTr$biozone.label=="K", 10) 

dim_height2<-PTr$dim_height+multiples # cumulative stratigraphic position of del13C on dimensionless time grid
PTr<-cbind(PTr, dim_height2)
PTr<-subset(PTr, dim_height2<12)


# =============================================================================
# Geographic and stratigraphic attributes
# =============================================================================


Iran<-subset(PTr, geographic.location == "Iran" )
China<-subset(PTr, geographic.location == "China" )

SB1<-5-(-303--1532)/(-303--2581) # sequence boundary

#--------------------------------------------------
# colors chronostratigraphic units of the 
# International Commision on Stratigraphy 
# (http://www.stratigraphy.org/index.php/ics-chart-timescale).
#--------------------------------------------------

triassic <- rgb(129,43,146, maxColorValue=255)
permian <- rgb(240,64,40, maxColorValue=255)
wuchiapingian <- rgb(252,180,162, maxColorValue=255)
changhsingian <- rgb(252,192,178, maxColorValue=255)
induan <- rgb(164,70,159, maxColorValue=255)


strat.colors.system <- c(triassic,permian) # system colors
strat.colors.stage <- c(induan, changhsingian, wuchiapingian) # stage colors

#---------------------------------------------------
# Chronostratigraphical bounds
#---------------------------------------------------

bounds.system <- c(7,11)
bounds.system.label <- c("Permian", "Triassic")
bounds.system.pos <- c(3.5, 9)


bounds.stage <- c(3, 7, 11)
bounds.stage.label <- c("Wuchiapingian", "Changhsingian", "Induan")
bounds.stage.pos <- c(1.5,5,9)


# =============================================================================
# Subsampling of del13C data by application of a sliding window
# =============================================================================



SampsizeBoot <- function(data,time, grid, windowsize){ 
  
  rarefrac_prep<-c()
  
  for(i in 1:(length(grid))){
    
    z<-length(data[time >= (grid[i]-(windowsize/2)) & time <=(grid[i]+(windowsize/2))])
    if(z==0){z<-NA}
    rarefrac_prep<-c(rarefrac_prep, z)
    
  }
  
  sampsize<-min(rarefrac_prep, na.rm=TRUE)
  return(sampsize)
} # determening the window with the lowest sample size. This value is subsequently used to determine the minimum sample size required for implementing the subsampling treatment.


#---------------------------------------------------
# Setting-up the grid for sliding time window subsampling
#---------------------------------------------------

grid.Iran    <- seq(0,10, 0.1)        # time grid Iran
grid.China   <- seq(2,10, 0.1)        # time grid China
windowsize   <- 1                     # size of sliding time window

#---------------------------------------------------

size.Iran<-SampsizeBoot(data=Iran$del13C,time=Iran$dim_height2, grid=grid.Iran, windowsize=windowsize) # Iran
size.China<-SampsizeBoot(data=China$del13C,time=China$dim_height2, grid=grid.China , windowsize=windowsize) # China

BootSlide <- function(n,data,time,grid,sampsize,windowsize,method){ 
  
  l0.boot <- matrix(NA, nrow=length(grid), ncol=n)
  
  for(i in 1:(length(grid))){
    
    f<-data[time >= (grid[i]-(windowsize/2)) & time <=(grid[i]+(windowsize/2))]
    
    if(length(f)==0){
      
      l0.boot[i,]<-rep(NA,n)
      
    } else {
      
      seR<-seq(0,1,0.025)
      
      for(t in 1:n){
        boot <- sample(f, sampsize, replace=TRUE)
        
        l0 <- method(data)
        if (is.na(l0)==TRUE) {
          m1   <- method(boot, na.rm=TRUE)}else{
            m1   <- method(boot)  
          }
        l0.boot[i,t] <- m1}      
    }  
  } 
  return(l0.boot)
}# this function represents the subsampling approach, giving both the complete set of subsamples as well as their corresponding summary statistics by adding the functions for the median, and the del13C value range of the interquartile range (IQR)-50% as well as the inter-percentile range (IPR)-95 % of the sample population.

#---------------------------------------------------
# Iran, subsampled and summary statistics
#---------------------------------------------------

BootMed.Iran<-BootSlide(n=999,data=Iran$del13C,time=Iran$dim_height2, grid=grid.Iran , windowsize=windowsize, sampsize=size.Iran, method=median)
BootIQR.Iran<-BootSlide(n=999,data=Iran$del13C,time=Iran$dim_height2, grid=grid.Iran , windowsize=windowsize, sampsize=size.Iran, method=function(x, seR=seq(0,1,0.25)){(quantile(x,seR))[4]-(quantile(x,seR))[2]})
BootIPR.Iran<-BootSlide(n=999,data=Iran$del13C,time=Iran$dim_height2, grid=grid.Iran , windowsize=windowsize, sampsize=size.Iran, method=function(x, seR=seq(0,1,0.025)){(quantile(x,seR))[40]-(quantile(x,seR))[2]})


BootSum.Iran<-cbind(rbind(BootMed.Iran, BootIQR.Iran, BootIPR.Iran), c(rep("Median", nrow(BootMed.Iran)), rep("IQR (50 %)", nrow(BootIQR.Iran)), rep("IQR (95%)", nrow(BootIPR.Iran) ))) # combined outcomes of summary statistics

#---------------------------------------------------
# China, subsampled and summary statistics
#---------------------------------------------------

BootMed.China<-BootSlide(n=999,data=China$del13C,time=China$dim_height2, grid=grid.China , windowsize=windowsize, sampsize=size.China, method=median)
BootIQR.China<-BootSlide(n=999,data=China$del13C,time=China$dim_height2, grid=grid.China , windowsize=windowsize, sampsize=size.China, method=function(x, seR=seq(0,1,0.25)){(quantile(x,seR))[4]-(quantile(x,seR))[2]})
BootIPR.China<-BootSlide(n=999,data=China$del13C,time=China$dim_height2, grid=grid.China , windowsize=windowsize, sampsize=size.China, method=function(x, seR=seq(0,1,0.025)){(quantile(x,seR))[40]-(quantile(x,seR))[2]})


BootSum.China<-cbind(rbind(BootMed.China, BootIQR.China, BootIPR.China), c(rep("Median", nrow(BootMed.China)), rep("IQR (50 %)", nrow(BootIQR.China)), rep("IQR (95%)", nrow(BootIPR.China) ))) # combined outcomes of summary statistics


# =============================================================================
# Construction of y-transect for visual weight (adapted from R code  reported in blogs by Felix Schonbrodt, http://www.nicebread.de/, and Solomon Hsiang, http://www.fight-entropy.com/). 
# Replication of the original license, 
# follows below;

# Copyright 2012 Felix Schonbrodt
# All rights reserved.
# 
# FreeBSD License
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#       
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#       
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER `AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# The views and conclusions contained in the software and documentation
# are those of the authors and should not be interpreted as representing
# official policies, either expressed or implied, of the copyright
# holder.

# Version history:
# 0.1: original code
# 0.1.1: changed license to FreeBSD; re-established compability to ggplot2 (new version 0.9.2)

## Visually weighted regression / Watercolor plots
## Idea: Solomon Hsiang, with additional ideas from many blog commenters

# =============================================================================


colorWeight<-function(data, grid, ylim=c(-10, 10), slices=1000)  {        
  
  # Compute median and CI limits of bootstrap
  
  CI.boot<- adply(data, 1, function(x) quantile(x, prob=c(.025, .5, .975, pnorm(c(-3, -2, -1, 0, 1, 2, 3))), na.rm=TRUE))[, -1]
  colnames(CI.boot)[1:10] <- c("LL", "M", "UL", paste0("SD", 1:7))
  CI.boot$x <- grid
  CI.boot$width <- CI.boot$UL - CI.boot$LL
  
  # Scale the CI width to the range 0 to 1 and flip it (bigger numbers = narrower CI)
  
  CI.boot$w2 <- (CI.boot$width - min(CI.boot$width))
  CI.boot$w3 <- 1-(CI.boot$w2/max(CI.boot$w2))
  
  # Convert bootstrapped spaghettis to long format
  
  b2 <- melt(data)
  b2$x <- grid
  colnames(b2) <- c("index", "B", "value", "x")
  
  
  
  # Range and number of tiles used in ggplot
  
  ylim <- c(-10, 10)
  slices <- 1000
  
  # Vertical cross-sectional density estimate
  
  d2 <- ddply(b2[, c("x", "value")], .(x), function(df) {
    res <- data.frame(density(df$value, na.rm=TRUE, n=slices, from=ylim[1], to=ylim[2])[c("x", "y")])
    
    colnames(res) <- c("y", "dens")
    return(res)
  }, .progress="text")
  
  maxdens <- max(d2$dens)
  mindens <- min(d2$dens)
  d2$dens.scaled <- (d2$dens - mindens)/maxdens   
  
  # Shading and coloring of density tiles
  
  shade.alpha=.1 # should the CI shading fade out at the edges? (by reducing alpha; 0 = no alpha decrease, 0.1 = medium alpha decrease, 0.5 = strong alpha decrease)
  
  
  # Tile approach
  
  d2$alpha.factor <- d2$dens.scaled^shade.alpha
  
  return(color=list(CI.boot, d2)) 
} # function attributing a value that weighs color to the confidence interval of the median trend lines constructed by multiple permutations of the subsampling routine from above.

#---------------------------------------------------
# Execution of the visual weight on the combined datasets 
#---------------------------------------------------

# Iran
colorMed.Iran<-colorWeight(data=BootMed.Iran, grid=grid.Iran)
colorIQR.Iran<-colorWeight(data=BootIQR.Iran, grid=grid.Iran)
colorIPR.Iran<-colorWeight(data=BootIPR.Iran, grid=grid.Iran)

CI.sum.Iran<-cbind(rbind((colorMed.Iran[[1]]), (colorIQR.Iran[[1]]), (colorIPR.Iran[[1]])), label=c(rep("Median", nrow((colorMed.Iran[[1]]))), rep("IQR (50 %)", nrow((colorIQR.Iran[[1]]))), rep("IPR (95%)", nrow((colorIPR.Iran[[1]])) )), stringsAsFactors=FALSE)
color.sum.Iran<-cbind(rbind((colorMed.Iran[[2]]), (colorIQR.Iran[[2]]), (colorIPR.Iran[[2]])), label=c(rep("Median", nrow((colorMed.Iran[[2]]))), rep("IQR (50 %)", nrow((colorIQR.Iran[[2]]))), rep("IPR (95%)", nrow((colorIPR.Iran[[2]])))), stringsAsFactors=FALSE)

# China
colorMed.China<-colorWeight(data=BootMed.China, grid=grid.China)
colorIQR.China<-colorWeight(data=BootIQR.China, grid=grid.China)
colorIPR.China<-colorWeight(data=BootIPR.China, grid=grid.China)

CI.sum.China<-cbind(rbind((colorMed.China[[1]]), (colorIQR.China[[1]]), (colorIPR.China[[1]])), label=c(rep("Median", nrow((colorMed.China[[1]]))), rep("IQR (50 %)", nrow((colorIQR.China[[1]]))), rep("IPR (95%)", nrow((colorIPR.China[[1]])) )), stringsAsFactors=FALSE)
color.sum.China<-cbind(rbind((colorMed.China[[2]]), (colorIQR.China[[2]]), (colorIPR.China[[2]])), label=c(rep("Median", nrow((colorMed.China[[2]]))), rep("IQR (50 %)", nrow((colorIQR.China[[2]]))), rep("IPR (95%)", nrow((colorIPR.China[[2]])))), stringsAsFactors=FALSE)



# =============================================================================
# Data plots (Figure 1)
# =============================================================================

theme=theme_set(theme_classic()) 

#---------------------------------------------------
# Stratigraphic del13C plot
#---------------------------------------------------

# Iran
I<-ggplot(colorMed.Iran[[1]], aes(x=x, y=M))+
  
  # Annotation for the extinction horizon
  
  annotate("text", label = "extinction \n horizon", x = 6, y = -5, size=3.5, color="black", angle=90)+
  geom_vline(xintercept=6, lty=2)+
  
  # Annotation for the Permian-Triassic boundary
  
  geom_vline(xintercept=7, lty=1)+
  
  # Annotation for the sea-level changes
  
  geom_rect(aes(xmax=SB1, xmin=0, ymax=10, ymin=8), fill="grey")+
  geom_rect(aes(xmax=11, xmin=SB1, ymax=10, ymin=8))+
  annotate("text", label = "lowstand", x = 2.5, y = 9, size=4, color="white")+
  annotate("text", label = "transgresssion", x = 7.5, y = 9, size=4, color="grey")+
  
  # Geographic location
  
  annotate("text", label= "Iran", x = 10, y = 7, size=4, fontface="bold")+
  
  # Data points
  
  geom_point(aes(x=dim_height2, y=del13C, colour=label), data=Iran, shape=1)+
  
  # Trend line  
  
  geom_line(color="black")+
  
  # Simulated seawater DIC-del13C line, requires times series simulation (CarbDiaModel.R)
  
  #geom_line(data=TimeSeries.sw.Iran[1:101,], aes(x=grid.Iran, y=d13C.sw), linetype=2, color="blue")+
  
  
  
  # Axis
  
  scale_y_continuous(limits=c(-10, 10), expand = c(0, 0))+
  ylab(expression(atop(paste(delta^13*C~"(\211 VPDB)"))))+
  labs(title= "(a)")+
  
  # Overwrites the labels for the aes x to create the names of the conodont zones 
  
  scale_x_discrete("biozones", limits=seq(0, 11,1) , expand = c(0, 0), labels=c("",paste(LETTERS[1:11])))+ 
  coord_cartesian(xlim=seq(0, 11,1))+
  
  theme(
    axis.line.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.line.y = element_line(color = "black"),
    legend.title = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    legend.text = element_text(size=5), 
    legend.position = c(0.25,0.26), 
    legend.key.size = unit(3,"mm"),
    legend.background = element_blank(),
    
    axis.text.x  = element_blank(), 
    axis.title.x  = element_blank(), 
    axis.text.y = element_text(size=7), 
    axis.title.y = element_text(size =7),
    plot.title = element_text(hjust=(0), size=9) )  




# China    
C<-ggplot(colorMed.China[[1]], aes(x=x, y=M))+
  
  # Annotation for the extinction horizon
  
  geom_vline(xintercept=6, lty=2)+
  
  # Annotation for the Permian-Triassic boundary
  
  geom_vline(xintercept=7, lty=1)+
  
  # Data points
  
  geom_point(aes(x=dim_height2, y=del13C, colour=label), data=China, shape=1)+
  
  # Trend line  
  
  geom_line(color="black")+
  
  # Simulated seawater DIC-d13C line, requires times series simulation (CarbDiaModel.R)
  
  #geom_line(data=TimeSeries.sw.China, aes(x=grid, y=d13C.sw), linetype=2, color="blue")+
  
  
  # Axis
  
  scale_y_continuous(limits=c(-10, 10), expand = c(0, 0))+
  ylab(expression(atop(paste(delta^13*C~"(\211 VPDB)"))))+
  labs(title= "(b)")+  
  
  # Geographic location
  
  annotate("text", label= "China", x = 10, y = 7, size=4, fontface="bold")+
  
  
  # Overwrites the labels for the aes x to create the names of the conodont zones 
  
  scale_x_discrete("biozones", limits=seq(0, 11,1) , expand = c(0, 0), labels=c("",paste(LETTERS[1:11])))+ 
  coord_cartesian(xlim=seq(0, 11,1))+
  
  
  # Chronological boxes for the systems
  
  geom_rect(aes(xmax=bounds.system[1], xmin=0, ymax=-9, ymin=-10), fill=strat.colors.system[2])+
  geom_rect(aes(xmax=bounds.system[2], xmin=bounds.system[1], ymax=-9, ymin=-10), fill=strat.colors.system[1])+
  
  annotate("text", label=bounds.system.label, x=bounds.system.pos, y=-9.5, size=3.5)+
  
  # Chronological boxes for the stages
  
  geom_rect(aes(xmax=bounds.stage[1], xmin=0,  ymax=-8, ymin=-9), fill=strat.colors.stage[3])+
  geom_rect(aes(xmax=bounds.stage[2], xmin=bounds.stage[1],  ymax=-8, ymin=-9), fill=strat.colors.stage[2])+
  geom_rect(aes(xmax=bounds.stage[3], xmin=bounds.stage[2],  ymax=-8, ymin=-9), fill=strat.colors.stage[1])+
  
  annotate("text", label=bounds.stage.label, x=bounds.stage.pos, y=-8.5, size=3.5)+
  
  theme(
    axis.line.x = element_line(color = "black"), 
    axis.line.y = element_line(color = "black"),
    legend.title = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    legend.text = element_text(size=5), 
    legend.position = c(0.245,0.33), 
    legend.key.size = unit(3,"mm"),
    legend.background = element_blank(),
    
    axis.text.x  = element_text(face = "bold", hjust = 5, size = 5), 
    axis.title.x  = element_text(size = 7), 
    axis.text.y = element_text(size=7), 
    axis.title.y = element_text(size =7),
    plot.title = element_text(hjust=(0), size=9))



# This guide transforms the aes-fill (sections) to two columns

I<-I+guides(colour=guide_legend(ncol=2))
C<-C+guides(colour=guide_legend(ncol=2))


# Code to override clipping

gI <- ggplot_gtable(ggplot_build(I))
gI$layout$clip[gI$layout$name == "panel"] <- "off"
grid.draw(gI)


gC <- ggplot_gtable(ggplot_build(C))
gC$layout$clip[gC$layout$name == "panel"] <- "off"
grid.draw(gC)


# Printing

grid.arrange(gI, gC, nrow=2) 

D13C<-arrangeGrob(gI, gC, nrow=2) 

ggsave("Figure1.pdf", D13C, height= 16, width=16, units="cm" )



# =============================================================================
# Visually weighted plots (Figures 3 and 4)
# =============================================================================


palette<-colorRampPalette(c("#a6bddb", "#034e7b"), bias=2)(20) # color palette

Ic<-ggplot(CI.sum.Iran, aes(x=x, y=M,  alpha=w3^3))+
  
  facet_grid(label~., scales="free")+
  
  # Annotation for the extinction horizon
  
  geom_vline(xintercept=6, lty=2)+
  
  # Annotation for the Permian-Triassic boundary
  
  geom_vline(xintercept=7, lty=1)+
  
  # Smoothed trend line
  
  geom_tile(data=color.sum.Iran, aes(x=x, y=y, fill=dens.scaled, alpha=alpha.factor))+
  scale_fill_gradientn("dens.scaled", colours=palette, guide=FALSE)+ 
  scale_alpha_continuous(range=c(0.001, 1), guide=FALSE)+
  geom_line(color="white")+
  
  # Axis  
  
  ylab(expression(atop(paste(delta^13*C~"(\211 VPDB)"))))+
  
  scale_y_continuous(limits=c(-10, 10), expand = c(0, 0))+
  
  # Overwrites the labels for the aes x to create the names of the conodont zones 
  
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



Cc<-ggplot(CI.sum.China, aes(x=x, y=M,  alpha=w3^3))+
  
  facet_grid(label~., scales="free")+
  
  # Annotation for the extinction horizon
  
  geom_vline(xintercept=6, lty=2)+
  
  # Annotation for the Permian-Triassic boundary
  
  geom_vline(xintercept=7, lty=1)+
  
  # Smoothed trend line
  
  geom_tile(data=color.sum.China, aes(x=x, y=y, fill=dens.scaled, alpha=alpha.factor))+
  scale_fill_gradientn("dens.scaled", colours=palette, guide=FALSE)+ 
  scale_alpha_continuous(range=c(0.001, 1), guide=FALSE)+
  geom_line(color="white")+
  
  # Axis  
  
  ylab(expression(atop(paste(delta^13*C~"(\211 VPDB)"))))+
  scale_y_continuous(limits=c(-10, 10), expand = c(0, 0))+
  
  # Overwrites the labels for the aes x to create the names of the conodont zones 
  
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


#---------------------------------------------------
# Biozone thickness and duration
#---------------------------------------------------





# Data frame conversion
med.THICK.Iran<-data.frame(median=apply(THICK[,c(1:7)],1, median, na.rm=TRUE), dim_height2=c(1:11), y=rep(1,11))
med.THICK.China<-data.frame(median=THICK[,8], dim_height2=c(1:11), y=rep(1,11))

palette=colorRampPalette(c("#54FF9F", "#2E8B57"), bias=2)(20) #color palette used in graph

# Iran thickness of biozones
tI<-ggplot(med.THICK.Iran, aes(x=dim_height2, y=y,  fill=median))+
  geom_tile()+
  scale_fill_gradientn("dens.scaled", colours=palette, guide=FALSE)+ 
  geom_rect(aes(xmin = 0.5, xmax = 11.5, ymin = 0.5, ymax = 1.5), fill = "transparent", color = "black", size = 0.5)+
  ylab("biozone \n thickness")+
  
  theme(
    axis.line.x = element_blank(), 
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size=5), 
    legend.position = c(0.18,0.29), 
    legend.key.size = unit(3,"mm"),
    legend.background = element_blank(),
    axis.title.y  = element_text(size=7), 
    axis.title.x  = element_blank(), 
    axis.text.x = element_blank(), 
    axis.text.y = element_blank())


# China thickness of biozones
tC<-ggplot(med.THICK.China, aes(x=dim_height2, y=y,  fill=median))+
  geom_tile()+
  scale_fill_gradientn("dens.scaled", colours=palette, guide=FALSE)+ 
  geom_rect(aes(xmin = 0.5, xmax = 11.5, ymin = 0.5, ymax = 1.5), fill = "transparent", color = "black", size = 0.5)+
  ylab("biozone \n thickness")+
  theme(
    axis.line.x = element_blank(), 
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size=5), 
    legend.position = c(0.18,0.29), 
    legend.key.size = unit(3,"mm"),
    legend.background = element_blank(),
    axis.title.y  = element_text(size=7), 
    axis.title.x  = element_blank(), 
    axis.text.x = element_blank(), 
    axis.text.y = element_blank())


# Data frame of biozones durations
DUR<-data.frame(DUR, y=rep(1, length(DUR)), x=c(1:length(DUR)))

D<-ggplot(DUR, aes(x=x, y=y,  fill=DUR))+
  geom_tile()+
  scale_fill_gradientn("dens.scaled", colours=palette, guide=FALSE)+ 
  geom_rect(aes(xmin = 0.5, xmax = 11.5, ymin = 0.5, ymax = 1.5), fill = "transparent", color = "black", size = 0.5)+
  ylab("biozone \n duration")+
  theme(
    axis.line.x = element_blank(), 
    axis.line.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size=5), 
    legend.position = c(0.18,0.29), 
    legend.key.size = unit(3,"mm"),
    legend.background = element_blank(),
    axis.text.x  = element_blank(), 
    axis.title.x  = element_blank(), 
    axis.text.y = element_blank(), 
    axis.title.y  = element_text(size=7))

empty <- ggplot()+geom_blank(aes(1,1))+
  theme(
    plot.background = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(), 
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

grid.arrange(arrangeGrob(arrangeGrob(empty, arrangeGrob(tI, D), empty, widths = c(1/100,98/100,1/100)), Ic, ncol=1, heights=c(1/6,5/6)),
             arrangeGrob(arrangeGrob(empty, arrangeGrob(tC, D), empty, widths = c(1/100,98/100,1/100)), Cc, ncol=1, heights=c(1/6,5/6)), ncol=2)

g<-arrangeGrob(arrangeGrob(arrangeGrob(empty, arrangeGrob(tI, D), empty, widths = c(1/25,23/25,1/25)), Ic, ncol=1, heights=c(1/6,5/6)),
               arrangeGrob(arrangeGrob(empty, arrangeGrob(tC, D), empty, widths = c(1/25,23/25,1/25)), Cc, ncol=1, heights=c(1/6,5/6)), ncol=2)



ggsave("Figure2.pdf", g,  height= 6, width=8, units="in")



# =============================================================================
# Comparative analyses of the (first-order) median trends recorded at 
# Meishan and Abadeh during multiple sampling expeditions (Figure 2)
# =============================================================================


C.int<-dlply(China, ~label , summarize, apply(BootSlide(n=999, data=del13C, time=dim_height2, grid=grid.China, windowsize=windowsize, method=median),1, FUN=median)) # obtaining trend lines for seperate studies on the Meishan profile

A.int<-dlply(Iran[which(Iran$locality=="Abadeh"),], ~label , summarize, apply(BootSlide(n=999, data=del13C, time=dim_height2, grid=grid.Iran, windowsize=windowsize, method=median),1, FUN=median)) # obtaining trend lines for seperate studies on the Abadeh profile


# Correlation of obtained firs-order trends Abadeh

k<-matrix(NA, nrow=length(unique((Iran[which(Iran$locality=="Abadeh"),])$label)), ncol=length(unique((Iran[which(Iran$locality=="Abadeh"),])$label)))

for(i in c(1:length(unique((Iran[which(Iran$locality=="Abadeh"),])$label)))){
  
  for(l in c(1:length(unique((Iran[which(Iran$locality=="Abadeh"),])$label)))){
    
    t<-data.frame(t1=unlist(A.int[[i]]), t2=unlist(A.int[[l]]))
    
    good<-complete.cases(t$t1,t$t2)
    
    if(nrow(t[good,])==0){ k[i,l]<-NA}else{
      
      k[i,l]<-summary(lm(t1~t2, data=t, singular.ok = TRUE))$r.squared }
    
    
  }}

k[k==0] <-NA 
k[k==1] <-NA

k<-data.frame(k,  row.names=unique((Iran[which(Iran$locality=="Abadeh"),])$label))
colnames(k)<-unique((Iran[which(Iran$locality=="Abadeh"),])$label)


k<-melt(k)
k<-cbind(k, name=rep(unique((Iran[which(Iran$locality=="Abadeh"),])$label), length(unique((Iran[which(Iran$locality=="Abadeh"),])$label))))


# Correlation of obtained firs-order trends Meishan

m<-matrix(NA, nrow=length(unique(China$label)), ncol=length(unique(China$label)))

for(i in c(1:length(unique(China$label)))){
  
  for(l in c(1:length(unique(China$label)))){
    
    t<-data.frame(t1=unlist(C.int[[i]]), t2=unlist(C.int[[l]]))
    
    good<-complete.cases(t$t1,t$t2)
    
    if(nrow(t[good,])==0){ m[i,l]<-NA}else{
      
      m[i,l]<-summary(lm(t1~t2, data=t, singular.ok = TRUE))$r.squared }
    
    
  }}

m[m==0] <-NA 
m[m==1] <-NA

m<-data.frame(m,  row.names=unique(China$label))
colnames(m)<-unique(China$label)

m<-melt(m)
m<-cbind(m, name=rep(unique(China$label), length(unique(China$label))))



# Comparative plot Abadeh 
Com.TA<-ggplot(k, aes(variable, name)) + geom_tile(aes(fill = value),colour = "white")+  
  # reshape plot to match size of Meishan plot
  geom_rect(aes(xmin=0.5,xmax=13.5, ymin=0.5, ymax=13.5), color="white", fill=NA)+ 
  # outline of plot area  
  geom_rect(aes(xmin=0.5,xmax=6.5, ymin=0.5, ymax=6.5), color="black", fill=NA)+
  scale_x_discrete(limits=unique(sort((Iran[which(Iran$locality=="Abadeh"),])$label)))+
  scale_y_discrete(limits=unique(sort((Iran[which(Iran$locality=="Abadeh"),])$label)))+
  scale_fill_gradient(expression("r"^2),low = "white",high = "steelblue", limits=c(0,1), breaks=c(0,0.5,1))+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=7),
        axis.text.y = element_text(size =7),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position=c(.7, .7),
        legend.key.width=unit(0.4, "line"),
        legend.key.height=unit(0.4, "line"),
        legend.text =element_text(size=7),
        axis.line = element_blank(),
        plot.margin=unit(c(0,0,1,0),"cm"))

# Comparative plot Meishan 
Com.TM<-ggplot(m, aes(variable, name)) + geom_tile(aes(fill = value),colour = "white")+  
  # outline of plot area 
  geom_rect(aes(xmin=0.5,xmax=13.5, ymin=0.5, ymax=13.5), color="black", fill=NA)+
  scale_x_discrete(limits=unique(sort(China$label)))+
  scale_y_discrete(limits=unique(sort(China$label)))+
  scale_fill_gradient(expression("r"^2),low = "white",high = "steelblue", limits=c(0,1), breaks=c(0,0.5,1))+
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size=7),
        axis.text.y = element_text(size =7),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",   
        axis.line = element_blank(),
        plot.margin=unit(c(0,0,1,0),"cm"))

# Equalizing plot heights

Abadeh<- ggplot_gtable(ggplot_build(Com.TA))  
Meishan<- ggplot_gtable(ggplot_build(Com.TM))  

Meishan$heights<- Abadeh$heights 

# Combining plots  


grid.arrange(Abadeh, Meishan, ncol=2)

# Printing plots

g<-arrangeGrob(Abadeh, Meishan, ncol=2)

ggsave("Figure3.pdf", g,  height= 8, width=18, units="cm")



####################################################################################################################################
# Rearranging data matrix (median trend lines) to fit input requirements for the time series 
# simulation by application of multiple reactive-transport model (CarbTrends.R)
####################################################################################################################################


ModelInput.Iran<-CI.sum.Iran[CI.sum.Iran$label=="Median",]

Age.unit<-outer((ModelInput.Iran$x>=0 & ModelInput.Iran$x<=1), (DUR[DUR$x==1,])$DUR)+
  outer((ModelInput.Iran$x>1 & ModelInput.Iran$x<=2), (DUR[DUR$x==2,])$DUR)+
  outer((ModelInput.Iran$x>2 & ModelInput.Iran$x<=3), (DUR[DUR$x==3,])$DUR)+
  outer((ModelInput.Iran$x>3 & ModelInput.Iran$x<=4), (DUR[DUR$x==4,])$DUR)+
  outer((ModelInput.Iran$x>4 & ModelInput.Iran$x<=5), (DUR[DUR$x==5,])$DUR)+
  outer((ModelInput.Iran$x>5 & ModelInput.Iran$x<=6), (DUR[DUR$x==6,])$DUR)+  
  outer((ModelInput.Iran$x>6 & ModelInput.Iran$x<=7), (DUR[DUR$x==7,])$DUR)+
  outer((ModelInput.Iran$x>7 & ModelInput.Iran$x<=8), (DUR[DUR$x==8,])$DUR)+
  outer((ModelInput.Iran$x>8 & ModelInput.Iran$x<=9), (DUR[DUR$x==9,])$DUR)+
  outer((ModelInput.Iran$x>9 & ModelInput.Iran$x<=10), (DUR[DUR$x==10,])$DUR)

L.unit<-outer((ModelInput.Iran$x>=0 & ModelInput.Iran$x<=1), (med.THICK.Iran[med.THICK.Iran$dim_height2==1,])$median)+
  outer((ModelInput.Iran$x>1 & ModelInput.Iran$x<=2), (med.THICK.Iran[med.THICK.Iran$dim_height2==2,])$median)+
  outer((ModelInput.Iran$x>2 & ModelInput.Iran$x<=3), (med.THICK.Iran[med.THICK.Iran$dim_height2==3,])$median)+
  outer((ModelInput.Iran$x>3 & ModelInput.Iran$x<=4), (med.THICK.Iran[med.THICK.Iran$dim_height2==4,])$median)+
  outer((ModelInput.Iran$x>4 & ModelInput.Iran$x<=5), (med.THICK.Iran[med.THICK.Iran$dim_height2==5,])$median)+
  outer((ModelInput.Iran$x>5 & ModelInput.Iran$x<=6), (med.THICK.Iran[med.THICK.Iran$dim_height2==6,])$median)+
  outer((ModelInput.Iran$x>6 & ModelInput.Iran$x<=7), (med.THICK.Iran[med.THICK.Iran$dim_height2==7,])$median)+
  outer((ModelInput.Iran$x>7 & ModelInput.Iran$x<=8), (med.THICK.Iran[med.THICK.Iran$dim_height2==8,])$median)+
  outer((ModelInput.Iran$x>8 & ModelInput.Iran$x<=9), (med.THICK.Iran[med.THICK.Iran$dim_height2==9,])$median)+
  outer((ModelInput.Iran$x>9 & ModelInput.Iran$x<=10), (med.THICK.Iran[med.THICK.Iran$dim_height2==10,])$median)

ModelInput.Iran<-cbind(ModelInput.Iran, IQR=CI.sum.Iran[CI.sum.Iran$label=="IQR (50 %)",], IPR=CI.sum.Iran[CI.sum.Iran$label=="IPR (95%)",],  Age.unit, L.unit )

ModelInput.China<-CI.sum.China[CI.sum.China$label=="Median",]

Age.unit<-outer((ModelInput.China$x>=0 & ModelInput.China$x<=1), (DUR[DUR$x==1,])$DUR)+
  outer((ModelInput.China$x>1 & ModelInput.China$x<=2), (DUR[DUR$x==2,])$DUR)+
  outer((ModelInput.China$x>2 & ModelInput.China$x<=3), (DUR[DUR$x==3,])$DUR)+
  outer((ModelInput.China$x>3 & ModelInput.China$x<=4), (DUR[DUR$x==4,])$DUR)+
  outer((ModelInput.China$x>4 & ModelInput.China$x<=5), (DUR[DUR$x==5,])$DUR)+
  outer((ModelInput.China$x>5 & ModelInput.China$x<=6), (DUR[DUR$x==6,])$DUR)+  
  outer((ModelInput.China$x>6 & ModelInput.China$x<=7), (DUR[DUR$x==7,])$DUR)+
  outer((ModelInput.China$x>7 & ModelInput.China$x<=8), (DUR[DUR$x==8,])$DUR)+
  outer((ModelInput.China$x>8 & ModelInput.China$x<=9), (DUR[DUR$x==9,])$DUR)+
  outer((ModelInput.China$x>9 & ModelInput.China$x<=10), (DUR[DUR$x==10,])$DUR)

L.unit<-
  
  outer((ModelInput.China$x>=2 & ModelInput.China$x<=3), (med.THICK.China[med.THICK.China$dim_height2==3,])$median)+
  outer((ModelInput.China$x>3 & ModelInput.China$x<=4), (med.THICK.China[med.THICK.China$dim_height2==4,])$median)+
  outer((ModelInput.China$x>4 & ModelInput.China$x<=5), (med.THICK.China[med.THICK.China$dim_height2==5,])$median)+
  outer((ModelInput.China$x>5 & ModelInput.China$x<=6), (med.THICK.China[med.THICK.China$dim_height2==6,])$median)+
  outer((ModelInput.China$x>6 & ModelInput.China$x<=7), (med.THICK.China[med.THICK.China$dim_height2==7,])$median)+
  outer((ModelInput.China$x>7 & ModelInput.China$x<=8), (med.THICK.China[med.THICK.China$dim_height2==8,])$median)+
  outer((ModelInput.China$x>8 & ModelInput.China$x<=9), (med.THICK.China[med.THICK.China$dim_height2==9,])$median)+
  outer((ModelInput.China$x>9 & ModelInput.China$x<=10), (med.THICK.China[med.THICK.China$dim_height2==10,])$median)

ModelInput.China<-cbind(ModelInput.China, IQR=CI.sum.China[CI.sum.China$label=="IQR (50 %)",], IPR=CI.sum.China[CI.sum.China$label=="IPR (95%)",],  Age.unit, L.unit)

####################################################################################################################################
####################################################################################################################################
