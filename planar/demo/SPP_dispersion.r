
## ----load, echo=FALSE,results='hide'-----------------------------------------
library(knitr)
library(ggplot2)
library(planar)
library(ggplot2)
require(reshape2)
library(gridExtra)
require(plyr)
opts_chunk$set(fig.path="sppdispersion/",
               warning=FALSE,error=FALSE,message=FALSE,tidy=FALSE)
library(ggplot2)
theme_set(theme_minimal() + theme(panel.border=element_rect(fill=NA)))


## ----kretschmann---------------------------------------------------------

k0 <- seq(1e-4, 3e-2, length=500)
wvl <- 2*pi/k0
silver <- epsAg(wvl)
gold <- epsAu(wvl)

dispersion <- function(material="silver"){
  material <- get(material)
  res <- recursive_fresnelcpp(k0=k0,
                              q=seq(0,1, length=500),
                              epsilon=list(1.5^2, material$epsilon, 1.0),
                              thickness=c(0, 50, 0),
                              polarisation='p')
  
  m <- melt(data.frame(k0=res$k0, R=res$R), id=c("k0"))
  m$q <- rep(Re(res$q), each=nrow(res$R))
  invisible(m)
}

m <- mdply(data.frame(material=c("silver","gold"), stringsAsFactors=FALSE), dispersion)

ggplot(m, aes(q, k0, fill=value)) +
  facet_wrap(~material, ncol=1) +
  geom_raster() + labs(fill = "R") +
  scale_x_continuous(expression(q==k[x] / k[1]), expand=c(0,0))+
  scale_y_continuous(expression(k[0]/nm^-1), expand=c(0,0))+
  theme_minimal()



## ----coupled-------------------------------------------------------------

k0 <- seq(1e-4, 2e-2, length=500)
wvl <- 2*pi/k0
silver <- epsAg(wvl)
gold <- epsAu(wvl)

coupled <- function(thickness = 50, material="silver"){
  material <- get(material)
  res <- recursive_fresnelcpp(k0=k0,
                              q=seq(1,1.4, length=500),
                              epsilon=list(1.0, material$epsilon, 1.0),
                              thickness=c(0, thickness, 0),
                              polarisation='p')
  
  m <- melt(data.frame(k0=res$k0, R=res$R), id=c("k0"))
  m$q <- rep(Re(res$q), each=nrow(res$R))
  invisible(m)
}

m <- mdply(data.frame(thickness=c(50, 1000)), coupled)

ggplot(m, aes(q, k0, fill=value)) +
  facet_wrap(~thickness, ncol=1) +
  geom_raster() + labs(fill = "R") +
  scale_x_continuous(expression(q==k[x] / k[1]), expand=c(0,0))+
  scale_y_continuous(expression(k[0]/nm^-1), expand=c(0,0))+
  theme_minimal()



