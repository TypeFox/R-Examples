## ----demo, message=FALSE, echo=FALSE-------------------------------------
library(knitr)

opts_chunk$set(cache=TRUE, cache.path="gaussian/",
               warning=FALSE,error=FALSE,message=FALSE,tidy=FALSE)

library(planar)
library(ggplot2)
library(grid)
library(plyr)


## ----pw, message=FALSE, echo=FALSE, fig.cap="Evolution of the near-field profile as a function of beam waist. The electric field is calculated at a fixed distance of 1nm from the metal surface. The (average) incident angle is kept constant, corresponding to the condition of optimum coupling for plane-wave illumination."----
struct <- list(wavelength=632.8,
               thickness=c(0,50,0),
               epsilon=c(1.5^2, epsAu(632.8)$epsilon, 1.0^2))

## first, check the plane wave result
pw <- multilayer(epsilon=as.list(struct$epsilon),
                 wavelength=struct$wavelength, thickness=struct$thickness, d=1,
                 angle=seq(0, pi/2, length=2e3), polarisation='p')

enhancement <- pw$Mr.perp[[2]] + pw$Mr.par[[2]]
maxi <- max(enhancement, na.rm=T)
spp <- pw$angle[which.max(enhancement)]


simulation <- function(w0=10){
  w0 <- w0*1e3
  xyz <- as.matrix(expand.grid(x=seq(-5*w0, 5*w0,length=100), y=0, z=51))
  
  res <- gaussian_near_field_ml(xyz, epsilon=struct$epsilon,
                                wavelength=struct$wavelength, 
                                thickness=struct$thickness,
                                w0=w0, alpha=spp, maxEval=1000)
  
  data.frame(xyz, field=res)
}

params <- data.frame(w0=c(10, 1e2, 1e3))
all <- mdply(params, simulation)

ggplot(all, aes(x/w0/1000, field, group=w0, colour=factor(w0)))+
  geom_line() + 
  geom_vline(aes(xintercept=0),lty=2) +
  geom_hline(aes(yintercept=maxi),lty=3) +
  annotate("text", label="plane-wave", y=maxi, x=-2.5, vjust=1, fontface="italic") +
  labs(x=expression(x/w[0]), y=expression("|E|"^2), 
       colour=expression(w[0]/mu*m)) +
  coord_cartesian(xlim=c(-5,5)) + theme_minimal()+
  guides(colour=guide_legend(reverse=TRUE)) +
  theme(panel.border=element_rect(fill=NA))


## ----weeber, echo=FALSE, fig.cap="Simulation of the electric field in the Kretschmann-Raether configuration. A gold layer (n=0.180 + 5.12i) surrounded by glass (n=1.52) and air is illuminated with TM-polarised light from the glass side. The electric field is calculated as a function of position along the interface (x) at a distance of 50nm (top panel) and 100nm (bottom panel) from the interface. The dotted curves correspond to the same simulation without a gold film."----

simulation <- function(d, probe=50, w0=10e3) {
  
  s <- list(wavelength=800, thickness=c(0,d,0),
                           epsilon=list(1.52^2, (0.180 + 5.12i)^2, 1.0^2))

  pw <- multilayer(epsilon=s$epsilon,
                   wavelength=s$wavelength, thickness=s$thickness, d=probe,
                   angle=seq(0, pi/2, length=2e3), polarisation='p')
  
  maxi <- max(pw$Mr.perp[[2]] + pw$Mr.par[[2]], na.rm=T)
  spp <- pw$angle[which.max(pw$Mr.perp[[2]] + pw$Mr.par[[2]])]
  
  xyz <- as.matrix(expand.grid(x=seq(-50e3, 150e3,length=100), y=0, z=s$thickness[2]+probe))
  
  res <- gaussian_near_field_layer(xyz, wavelength=s$wavelength,
                                epsilon=unlist(s$epsilon), thickness=s$thickness,
                                w0=w0, alpha=spp, maxEval=5000)
  data.frame(xyz, field=res/max(res))
}

params <- data.frame(d=c(50, 100))
all <- mdply(params, simulation)

bare <- simulation(0, 50) # no metal layer

peak <- function(d){
  peak <- subset(d, field == max(field))
  peakx <- peak$x /1000
  peakl <- round(peakx,2)
  peaky <- peak$field
  data.frame(peakx=peakx, peakl=peakl, peaky=peaky)
}

ann <- ddply(all, "d", peak)

ggplot(all, aes(x/1000, field))+ facet_grid(d~., scales="free")+
  geom_line()  +
  geom_line(data=bare, linetype="dotted") +
  geom_vline(aes(xintercept=0),lty=2) +
  geom_blank(data=ann, aes(y=peaky*1.1, x=peakx)) +
  geom_text(data=ann, aes(label=peakl, y=peaky, x=peakx), hjust=0, vjust=0, fontface="italic") +
  scale_y_continuous(expand=c(0,0)) +
  labs(x=expression(x/mu*m), y=expression("(normalised) |E|"^2)) +
  guides(colour="none") +
  theme_minimal()+
  theme(panel.border=element_rect(fill=NA)) +
  theme()

