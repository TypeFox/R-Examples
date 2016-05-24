## ----load,message=FALSE--------------------------------------------------
library(cda)
library(rgl)
library(ggplot2)
library(reshape2)
library(plyr)

## ----cluster, rgl=TRUE,echo=-12,tidy=FALSE,fig.width=3,fig.height=3,fig.path="basic-"----

# dielectric function
wvl <- seq(400, 900)
gold <- epsAu(wvl)

# define a cluster of particles
cl <- list(r = rbind(c(0, 0, 0),
                      c(0, 0, -100)),
            angles = rbind(c(0, 0, 0),
                           c(pi/4, 0, 0)),
            sizes = rbind(c(30, 10, 10),
                          c(30, 10, 10)))

# visualise
rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="gold")
rgl.viewpoint( theta = 0, phi = 20, fov = 70, zoom = 1)
rgl_annotate()



## ----linear,echo=TRUE,tidy=FALSE,fig.path="basic-", fig.height=4---------

linear <- dispersion_spectrum(cl, gold)
ggplot(linear, aes(wavelength, value, linetype=type)) +
  facet_wrap(~polarisation) + geom_path()


## ----oa,echo=TRUE,tidy=FALSE,fig.path="basic-",fig.width=8---------------
circular <- circular_dichroism_spectrum(cl, gold)

ggplot(circular, aes(wavelength, value, color=variable)) + 
  facet_grid(type~., scales="free") + geom_line()



