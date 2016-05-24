
## ----load,message=FALSE--------------------------------------------------
library(cda)
library(rgl)
library(ggplot2)
library(reshape2)
library(plyr)

## ----cluster, rgl=TRUE,echo=-12,tidy=FALSE,fig.width=3,fig.height=3,fig.path="dimercd-"----

# dielectric function
wvl <- seq(400, 900)
gold <- epsAu(wvl)

# define a dimer with dihedral angle
cl <- cluster_dimer(d=100, 
                    dihedral=45*pi/180, alpha1=10*pi/180, alpha2=0,
                    a=35, b=12)

# visualise
rgl.ellipsoids(cl$r, cl$sizes, cl$angles, col="gold")
rgl.viewpoint( theta = 0, phi = 20, fov = 70, zoom = 1)
rgl_annotate()



## ----comparison,echo=TRUE,tidy=FALSE,fig.path="dimercd-",fig.width=8-----
  
dimer <- function(dihedral=45, ...){
  cl <- cluster_dimer(dihedral = dihedral * pi/180, ...)
  circular_dichroism_spectrum(cl, material = gold)
  
}
params <- data.frame(dihedral=c(45, 30, 10, 0, -10, -30,-45),
                     d=100, a=35, b=12, 
                     alpha1=10*pi/180, alpha2=0)
comparison <- mdply(params, dimer)

p <- 
  ggplot(data=comparison) + 
  facet_grid(type~variable, scales="free") +
  geom_line(aes(wavelength, value, 
                colour=factor(dihedral))) +
  labs(y=expression(sigma*" /"*nm^2),
       x=expression(wavelength*" /"*nm), colour="dihedral angle") +
         scale_colour_brewer(type="div", palette=3)

p


