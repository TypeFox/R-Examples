
## ----load,message=FALSE--------------------------------------------------
library(cda)
library(rgl)
library(ggplot2)
library(reshape2)
library(plyr)

## ----cluster, rgl=TRUE,echo=-9,tidy=FALSE,fig.width=3,fig.height=3,fig.path="chain-"----
# dielectric function
wvl <- seq(400, 900)
gold <- epsAu(wvl)

cl <- cluster_chain(N=10, pitch=500)
rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="gold")
# visualise
rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="gold")
rgl.viewpoint( theta = 0, phi = 20, fov = 70, zoom = 1)


## ----comparison,echo=TRUE,tidy=FALSE,fig.path="chain-"-------------------

chain <- function(N, pitch = 500, ...){
  cl <- cluster_chain(N, pitch, ...)
  dispersion_spectrum(cluster = cl, material = gold)
}
  
params <- data.frame(N=c(1, 10, 50))
comparison <- mdply(params, chain)

p <- 
  ggplot(data=comparison)+ facet_wrap(~type, ncol=1, scales="free")+
labs(y=expression(sigma*" /"*nm^2),
       x=expression(wavelength*" /"*nm),
       colour = expression(N), linetype=expression(polarisation))+
  geom_line(aes(wavelength, value/N, linetype=polarisation,
                colour=factor(N),
                group=interaction(N,polarisation)))

p


