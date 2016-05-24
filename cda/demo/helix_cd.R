
## ----load,message=FALSE--------------------------------------------------
library(cda)
library(rgl)
library(ggplot2)
library(reshape2)
library(plyr)

## ----cluster, rgl=TRUE,echo=-12,tidy=FALSE,fig.width=3,fig.height=3,fig.path="helix-"----

# dielectric function
wvl <- seq(400, 900)
gold <- epsAu(wvl)

# define a helix
cl1 <- cluster_helix(N=20, R0=500, pitch=1000, 
                          delta=pi/7, delta0=0, right=TRUE,
                          a=100, b=50, c=50,
                          angles='helix')
cl2 <- cluster_helix(N=20, R0=500, pitch=1000, 
                          delta=pi/7, delta0=0, right=TRUE,
                          a=50, b=50, c=50,
                          angles="helix")
hel <- helix(N = 20, R0 = 500, pitch = 1000, delta = pi/7, 
             delta0 = 0, right = TRUE)
# visualise
rgl.ellipsoids(cl1$r, cl1$sizes, cl1$angles, col="gold")
lines3d(hel$smooth, lwd=1, col="red")
  shift <- cbind(rep(1000, nrow(cl2$r)),0,0)
  shifts <- cbind(rep(1000, nrow(hel$smooth)),0,0)
rgl.ellipsoids(cl2$r+shift, cl2$sizes, cl2$angles, col="gold")
lines3d(hel$smooth + shifts, lwd=1, col="red")


## ----comparison,echo=TRUE,tidy=FALSE,fig.path="helix-",fig.width=8-------
  
simulation <- function(N=3, ar=1, ...){
  cl <- cluster_helix(N, R0=12, pitch=15, 
                          delta=pi/2, delta0=0, right=TRUE,
                          a=5, b=5/ar, c=5/ar,
                          angles="helix")
  circular_dichroism_spectrum(cl, material = gold, medium=1.33)
  
}
params <- expand.grid(N=seq(2, 7), ar= c(1, 1.1))
comparison <- mdply(params, simulation)

p <- 
  ggplot(data=subset(comparison, type == "CD" & variable == "extinction")) + 
  facet_grid(ar ~ ., scales="free") +
  geom_line(aes(wavelength, value/N, 
                colour=factor(N))) +
  labs(y=expression(sigma[ext]*" /"*nm^2),
       x=expression(wavelength*" /"*nm), colour="# particles") 

p


