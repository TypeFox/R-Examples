# geosptDemo
fanfare <- function(stuff, cexx=2.5) {
   plot(c(0,1),c(0,1),type="n",axes=FALSE, xlab="", ylab="")
   text(0.5,0.5, stuff, cex=cexx)
 }

old.par <- par(no.readonly = TRUE)

par(ask=TRUE)

fanfare("I. Geostatistical modelling")
# graph.rbf function
data(preci)
coordinates(preci)<-~x+y
# optimizing eta
graph.rbf(prec~1, preci, eta.opt=TRUE, rho.opt=FALSE, n.neigh=9, func="TPS", np=40, eta.dmax=0.2, P.T=TRUE)

# pocket.plot function
library(gstat) 
data(coalash) 
plot(coalash[,1:2], type="n", xlab="x", ylab="y") 
text(coalash$x,coalash$y,coalash$coalash,cex=0.6)

# Pocket plot in the north-south direction. 
# Units on the vertical axis are root (\% coal ash) 

# Plot generated with the function pocket.plot 
# Clearly rows 2, 6, and 8 are atypical 

# This serves as verification that these rows are potentially problematic

# Analysis of local stationarity in probabilities of the coal in south-north direction 

pocket.plot(coalash, "PPR", coalash$x, coalash$y, coalash$coalash, FALSE)


fanfare("II. Design of Optimal Spatial\n Sampling Networks", 2.3)

# simPtsOptNet function

data(COSha30)
data(COSha30map)
data(lalib)

## Calculate the sample variogram for data, generate the variogram model and  
## fit ranges and/or sills from the variogram model to the sample variogram
ve <- variogram(CorT~1, loc=~x+y, data=COSha30, width = 236.0536)
PSI <- 0.0001531892; RAN <- 1031.8884; NUG <- 0.0001471817
m.esf <- vgm(PSI, "Sph", RAN, NUG)
(m.f.esf <- fit.variogram(ve, m.esf))

## Number of additional points to be added to the network
npoints <- 5

## Optimize the location of the additional points
## Users can visualize how the location of the additional points is optimized if plotMap is set to TRUE
plot(lalib)
par(ask=FALSE)
optnets <- simPtsOptNet(CorT~ 1, loc=~ x+y, COSha30, m.f.esf, n=npoints, 
    popSize=30, generations=25, xmin=bbox(lalib)[1], ymin=bbox(lalib)[2], 
    xmax=bbox(lalib)[3], ymax=bbox(lalib)[4], plotMap=TRUE, spMap=lalib)

par(ask=TRUE)
## Summary of the genetic algorithm results
summary(optnets, echo=TRUE)

## Graph of best and mean evaluation value of the objective function 
## (average standard error) along generations
plot(optnets)

## Find and plot the best set of additional points (best chromosome) in   
## the population in the last generation
(bnet <- bestnet(optnets))
l1 = list("sp.polygons", lalib)
l2 = list("sp.points", bnet, col="green", pch="*", cex=5)
spplot(COSha30map, "var1.pred", main="Location of the optimized points", 
    col.regions=bpy.colors(100), scales = list(draw =TRUE), xlab ="East (m)", 
    ylab = "North (m)", sp.layout=list(l1,l2))

### Average standard error of the optimized additional points
min(optnets$evaluations)
par(old.par)


