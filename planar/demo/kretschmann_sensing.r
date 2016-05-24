
## ----load, echo=FALSE,results='hide'-----------------------------------------
library(knitr)
library(ggplot2)
opts_chunk$set(fig.path="kretschmannsensing/",
               warning=FALSE,error=FALSE,message=FALSE,tidy=TRUE)
library(ggplot2)
theme_set(theme_minimal() + theme(panel.border=element_rect(fill=NA)))


## ----setup, results='hide'----------------------------------------------------
library(planar)
library(ggplot2)
require(reshape2)
library(gridExtra)
require(plyr)



## ----simulation----------------------------------------------------------
wvl <- 632.8
gold <- epsAu(wvl)
simulation <- function(thickness = 0, n1 = n2, n2 = 1.0){
  results <- recursive_fresnelcpp(epsilon=list(1.5^2, gold$epsilon, n1^2, n2^2),
                                  wavelength=gold$wavelength, 
                                  thickness=c(0, 50, thickness, 0),
                                  angle=seq(0, pi/2, length=2e3), polarisation='p')
  data.frame(results[c("angle", "R")])
}


## ----loop, fig.width=12--------------------------------------------------

## loop over parameters
parameters <- function(res=10) 
  data.frame(n2 = seq(1.0, 1.5, length=res))

d1 <- mdply(parameters(10), simulation)
d2 <- mdply(parameters(300), simulation)

p1 <- 
ggplot(d1) +
  geom_line(aes(angle*180/pi, R, colour=n2, group=n2)) +
  scale_y_continuous("Reflectivity", expand=c(0,0), limits=c(0,1))+
  scale_x_continuous("Internal angle /degrees", expand=c(0,0), 
                     breaks=seq(0,90, by=15)) +
  guides(colour=guide_legend())

## colour map
p2 <- 
ggplot(d2) +
  geom_raster(aes(angle*180/pi, n2, fill=R)) +
  scale_y_continuous("n", expand=c(0,0))+
  scale_x_continuous("Internal angle /degrees", expand=c(0,0), 
                     breaks=seq(0,90, by=15)) 

grid.arrange(p1, p2, nrow=1)


## ----variation-----------------------------------------------------------

## loop over parameters
parameters <- function(res=10) 
  expand.grid(thickness = c(10, 50, 100, 400),
              n1 = seq(1.0, 1.5, length=res))

d1 <- mdply(parameters(10), simulation, n2=1.0)

ggplot(d1) + facet_grid(thickness ~ . , scales = "free") + 
  geom_line(aes(angle*180/pi, R, colour=n1, group=n1)) +
  scale_y_continuous("Reflectivity", expand=c(0,0), limits=c(0,1))+
  scale_x_continuous("Internal angle /degrees", expand=c(0,0), 
                     breaks=seq(0,90, by=15)) +
  guides(colour=guide_legend()) +
  theme(panel.margin = unit(1, "lines"))


