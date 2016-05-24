
## ----load, echo=FALSE,results='hide'-----------------------------------------
library(knitr)
library(ggplot2)
opts_chunk$set(fig.path="kretschmannreflectivity/", fig.width=10,
               warning=FALSE,error=FALSE,message=FALSE,tidy=TRUE)
library(ggplot2)
theme_set(theme_minimal() + theme(panel.border=element_rect(fill=NA)))
library(RColorBrewer)
col <- brewer.pal(3,"PRGn")


## ----setup, results='hide'----------------------------------------------------
library(planar)
library(ggplot2)
require(reshape2)
library(gridExtra)
require(plyr)


## ----simulation------------------------------------------------------------------------
wvl <- 632.8
gold <- epsAu(wvl)
results <- recursive_fresnelcpp(epsilon=list(1.5^2, gold$epsilon, 1.0),
                                wavelength=gold$wavelength, thickness=c(0, 50, 0),
                                angle=seq(0, pi/2, length=2e3), polarisation='p')
str(results)


## ----reflectivity,fig.width=10-------------------------------------------
m <- data.frame(results[c("angle", "R")])

tir <- asin(1/1.5) * 180/pi
 
ggplot(m) +
  geom_vline(aes(xintercept=x),
             data=data.frame(x=tir),
             linetype=2,color="grey50") +
  geom_line(aes(angle*180/pi, R)) +
  scale_y_continuous("Reflectivity", expand=c(0,0), limits=c(0,1))+
  scale_x_continuous("Internal angle /degrees", expand=c(0,0), 
                     breaks=seq(0,90, by=15)) 
  


## ----loop----------------------------------------------------------------

simulation <- function(thickness = 50){
results <- recursive_fresnelcpp(epsilon=list(1.5^2, gold$epsilon, 1.0^2),
                                wavelength=gold$wavelength, 
                                thickness=c(0, thickness, 0),
                                angle=pi/180*seq(15, 60, length=500), 
                                polarisation='p')
data.frame(results[c("angle", "R")])

}

## loop over parameters
parameters <- function(res=10) 
  data.frame(thickness = seq(0, 100, length=res))

d1 <- mdply(parameters(10), simulation)
d2 <- mdply(parameters(300), simulation)


p1 <- 
ggplot(d1) +
  geom_line(aes(angle*180/pi, R, colour=thickness, group=thickness)) +
  scale_y_continuous("Reflectivity", expand=c(0,0), limits=c(0,1))+
  scale_x_continuous("Internal angle /degrees", expand=c(0,0), 
                     breaks=seq(0,90, by=15)) +
  guides(colour=guide_legend()) 

## colour map
p2 <- 
ggplot(d2) +
  geom_raster(aes(angle*180/pi, thickness, fill=R)) +
  scale_y_continuous("thickness", expand=c(0,0))+
  scale_x_continuous("Internal angle /degrees", expand=c(0,0), 
                     breaks=seq(0,90, by=15)) 

grid.arrange(p1, p2, nrow=2)


## ----variation-----------------------------------------------------------
minimum <- ddply(subset(d2, angle > tir*pi/180 & thickness > 5), .(thickness), summarize, 
                 angle = angle[which.min(R)] * 180/pi,
                 min = min(R))
ggplot(melt(minimum, id="thickness")) + 
  facet_grid(variable~., scales="free") +
  geom_line(aes(thickness, value)) +
  labs(y="", x="thickness /nm")
  


