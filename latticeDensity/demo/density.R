require(latticeDensity)


## LatticeDensity is designed to produce density maps
## accountaing for the effects of boundaries and holes.
## The minimum data needed to run latticeDensity is two
## 2-column matrices.  One gives the vertices of a polygon
## and the other the Easting and Northing locations of the
## observations.  Note that an observation outside the boundary
## polygon is moved automatically to the nearest location 
## inside the polygon.  



plot.new()
data(polygon1)
nodeFillingOutput = nodeFilling(poly=polygon1,node.spacing=0.02)
plot(nodeFillingOutput)

## The function nodeFillingOutput fills the polygon with 
## a rectangular pattern of nodes, here at spacing 0.01.
## We should examine the plot to determine whether the
## density is high enough to fill the polygon and limn
## interesting parts of the boundary.  If we have a list
## of polygonal holes, it can be used as an arugment in
## nodeFilling.

formLatticeOutput = formLattice(nodeFillingOutput)
plot(formLatticeOutput)

##  formLatticeOutput connects neighboring nodes into a
##  lattice.  Note to see whether there are connections
##  where there should not be (for instance, across a
##  causeway) or whether some connections are missing.
##  You can edit the neighbor structure using the
##  function editLattice if need be.

Pointdata = csr(polygon1,100)
Pointdata = Pointdata[Pointdata[,1]<0.5,]
plot(rbind(polygon1,polygon1[1,]),type="l")
points(Pointdata,pch=19)

out = crossvalDensity(formLatticeOutput,PointPattern=Pointdata, 
  M=0.5,num.steps = 150)

##  crossvalDensity uses crossvalidation to choose an optimal
##  path length k.  The larger k, the smoother the map.
##  M is the probability that the random walk remains in
##  place at each step and thus also governs the smoothness
##  of the density map.  The function also plots the ucv
##  criterion vs number of steps.  Note that in some simulations
##  ones gets a plot of ucv vs k in which the ucv is still decreasing
##  at the right side of the plot.  In this case we could rerun 
##  crossvalDensity with a larger num.steps.  This simulation is prone
##  to have a ucv vs k curve that is very flat for most values of k,
##  since the distribution of observations is homogenous Poisson to
##  the left of the causeway.

densityOut = createDensity(formLatticeOutput,
  PointPattern=Pointdata, k=out[[2]],intensity=FALSE, sparse = TRUE)
plot(densityOut)

##  Once the optimal number of steps in the diffusion is selected (either
##  by eye or by the crossvalDensity function), the function createDensity
##  produces an object suitable for constructing a contour plot.

homerange(densityOut, percent = 0.95)



## The homerange function fills in the minimal area in the
## region that has integrated density of 0.95 or higher.
devAskNewPage(ask = FALSE)
