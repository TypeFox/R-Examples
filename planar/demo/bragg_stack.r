
## ----setup, echo=FALSE,results='hide'------------------------------------
library(ggplot2)
library(knitr)
opts_chunk$set(fig.path="braggstack/",
               warning=FALSE,error=FALSE,message=FALSE, tidy=FALSE)
theme_set(theme_minimal() + theme(panel.border=element_rect(fill=NA)))


## ----stack, results='hide'-----------------------------------------------
library(planar)

make_stack <- function(n = 3, wavelength=seq(300, 700),
                      lambda0 = 500, thickness = lambda0/4, 
                      nH=2.3, nL=1.38, nS=1.52,
                      angle=0){

  epsilon.list <- as.list(c(1, rep(c(nL, nH)^2, n), nS))
  thickness.list <- c(0, rep(thickness/c(nL, nH), n), 0)
  
  res <- recursive_fresnelcpp(wavelength=wavelength, angle=angle,
                       epsilon=epsilon.list, thickness=thickness.list,
                       polarisation='p')
data.frame(res[c("wavelength","k0", "R")])
}


## ----simulation----------------------------------------------------------
require(plyr)

params <- expand.grid(n=c(2, 4, 8), angle=seq(0,pi/3, by=pi/6))
all <- mdply(params, make_stack)
all$angle.d <- 180/pi*all$angle

library(ggplot2)
ggplot(all)+ facet_grid(angle.d ~ .)+
  geom_line(aes(2*pi/k0, R, colour=factor(n))) +
  labs(colour="layers") +
  scale_x_continuous("wavelength /nm")+
  scale_y_continuous("Reflectivity", expand=c(0,0), lim=c(0,1))



