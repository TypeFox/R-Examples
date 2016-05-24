
## ----load,message=FALSE, echo=1:6----------------------------------------
require(cda)
require(ggplot2)
library(reshape2)
library(plyr)

## ----cluster, tidy=FALSE, fig.path='averaging-'--------------------------
gold <- epsAu(seq(600, 800))

cl <- cluster_dimer(d=100, 
              dihedral=10*pi/180, alpha1=20*pi/180, alpha2=0,
              a=35, b=12)

# achiral cluster (plane of symmetry)
cl2 <- cluster_dimer(d=100, 
              dihedral=0*pi/180, alpha1=20*pi/180, alpha2=0,
              a=35, b=12)


## ----comparison, tidy=FALSE, fig.path='averaging-'-----------------------
params <- expand.grid(Nquad=c(10, 50, 100, 1000),
                       averaging=c("grid", "GL", "QMC"),
                       stringsAsFactors=FALSE)

comparison <- mdply(params, circular_dichroism_spectrum, cluster=cl, material=gold)
cheap <- circular_dichroism_spectrum(cluster=cl, material=gold, averaging="cheap")
converged <- circular_dichroism_spectrum(cluster=cl, material=gold, averaging="QMC", Nquad=5000)

p <- 
  ggplot(subset(comparison, type == "CD" & variable == "extinction"),
         aes(wavelength, value)) + 
  facet_grid(averaging~., scales="free")+
  geom_path(aes(colour=factor(Nquad), group=Nquad))+
  geom_path(data=subset(cheap, type == "CD" & variable == "extinction"), linetype=2)+
  geom_path(data=subset(converged, type == "CD" & variable == "extinction"))+
  labs(y=expression(sigma[ext]/nm^2),
       x=expression(wavelength/nm), colour=expression(N))

p


## ----achiral, tidy=FALSE, fig.path='averaging-'--------------------------
params <- expand.grid(Nquad=c(10, 100, 1000, 5000),
                       averaging=c("grid", "GL", "QMC", "cheap"),
                       stringsAsFactors=FALSE)

comparison <- mdply(params, circular_dichroism_spectrum, cluster=cl2, material=gold)

p <- 
  ggplot(subset(comparison, type == "CD" & variable == "extinction")) + 
  facet_grid(averaging~.)+
  geom_path(aes(wavelength, value, colour=factor(Nquad), group=Nquad))+
  labs(y=expression(sigma[ext]/nm^2),
       x=expression(wavelength/nm), colour=expression(N))

p


