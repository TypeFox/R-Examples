
## ----load, echo=FALSE,results='hide', message=FALSE--------------------------
library(knitr)
library(planar)
library(ggplot2)
require(reshape2)
require(plyr)
opts_chunk$set(fig.path="decayrates/", 
               warning=FALSE,error=FALSE,message=FALSE,tidy=TRUE)
library(ggplot2)
theme_set(theme_minimal() + theme(panel.border=element_rect(fill=NA)))


## ----setup, results='hide'----------------------------------------------------

wvl <- seq(200, 1000,by=2)
silver <- epsAg(wvl)
gold <- epsAu(wvl)

distance <- function(d=1, material="silver", ...){

  material <- get(material)
  
  dl <- dipole(d=d,
               wavelength = material$wavelength,
               epsilon = list(incident=1.0^2, material$epsilon),
               thickness = c(0, 0),
               Nquadrature1 = 1e3, Nquadrature2 = 5e3, GL = FALSE,
               Nquadrature3 = 5e3, qcut = NULL, rel.err=1e-3,
               show.messages = FALSE)

  message(attr(dl, "comment"))
  
  m <- melt(dl, id = "wavelength")

  m$orientation <- m$variable
  
  levels(m$orientation) <- list(perpendicular="Mtot.perp",
                                perpendicular="Mrad.perp",
                                parallel="Mtot.par",
                                parallel="Mrad.par")
  
  levels(m$variable) <- list(Mtot="Mtot.perp",
                             Mtot="Mtot.par",
                             Mrad="Mrad.perp",
                             Mrad="Mrad.par")
  invisible(m)
}


## ----simulation, message=TRUE, fig.width=10------------------------------
params <- expand.grid(d=c(1,5,10), material=c("silver", "gold"), 
                      stringsAsFactors = FALSE)
all <- mdply(params, distance)

ggplot(all, aes(wavelength, value, colour=factor(d), linetype=orientation))+
  facet_grid(variable~material, scales="free_y") + 
  geom_path() + labs(colour="distance /nm", y="EM enhancement factor", x="wavelength /nm")+
  scale_y_log10() 


