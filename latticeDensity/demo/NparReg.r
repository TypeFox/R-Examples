require(latticeDensity)

#  The functions addQuantVar, createNparReg and crossvalNparReg
#  use the kernels based on random walks on lattices to make
#  kernel non-parametric regression maps in irregular regions.
#
#  The first step in this demo is opening the dataset nparExample, which
#  constains a boundary polygon called polygon2 and a grid of response
#  locations called grid2.  Note that the trickiest part of running functions
#  in the latticeDensity package is getting the polygons correct.  The vertices
#  have to be in order (either clockwise or counterclockwise).  Many polygons
#  such as extracted from shapefiles are defined in pieces, so that the
#  points, when plotted in the order found in the matrix, do not form the
#  correct polygon.

data(nparExample)
attach(nparExample)
plot.new()


#  For this dataset, we will simulate a response variable, with a value
#  simulated at each locations in grid2.  Note that the responses do not
#  have to be measured on a grid, any set of locations withing the boundary
#  will do.  In fact, even locations outside the boundary (due to measurement
#  error, perhaps) will be mapped to the nearest locations within the boundary.

index1 = (grid2[,2]<0.8)|(grid2[,1]>0.6)
Z = rep(NA,length(grid2[,1]))
n1 = sum(index1)
n2 = sum(!index1)
Z[index1] = 3*grid2[index1,1] + 4 + rnorm(n1,0,sd=0.4)
Z[!index1] = -2*grid2[!index1,1] + 4 + rnorm(n2,0,sd=0.4)

#  Next, I'll plot the values of the response variable at the appropriate
#  locations in the polygon.  Note the trends in the values of the response.

plot(rbind(polygon2,polygon2[1,]),type="l")
points(grid2,pch=19,cex=0.5,xlim=c(-0.1,1))
text(grid2,labels=round(Z,1),pos=4,cex=0.5)

#  The nodeFilingOutput function fills the region with nodes that are spaced
#  according (unsurprisingly) the argument node.spacing.  You want the nodes
#  sufficiently close together so that measured responses don't shift too much
#  when moved to the nearest node (the red and green dots show the response
#  and the node it is shifted to) and all the major features of the polygon
#  are filled with nodes.  The only upper limit is computing time and memory.

nodeFillingOutput = nodeFilling(poly=polygon2,node.spacing=0.025)
plot(nodeFillingOutput)

#  formLatticeOutput connects the nodes into a lattice upon which the 
#  random walk moves.  You want to make sure that you don't have (1) nodes
#  that should not be connected (i.e. on either side of a causeway) with
#  a link, and (2) nodes that should be connected lack a link.  If the 
#  plot looks fine, you can go straight to crossvalNparReg or createNparReg.
#  Otherwise, you can edit the lattice using the function editLattice.
#  You may have to magnify your screen, as editLattice uses mouse clicks to
#  add or remove links.

formLatticeOutput = formLattice(nodeFillingOutput)
plot(formLatticeOutput)

#  Next we use the function crossvalNparReg to use least-squares
#  crossvalidation to pick the optimal number of steps (smoothing).
#  If you are brave, you can skip crossvalidation and select a value
#  of the smoothing parameter k that seems to give good plots.  The
#  response variable is in Z, the locations of the responses are in 
#  PointPattern, M is another smoothing parameter (usually left at 0.5),
#  num.steps is the maximum number of steps examined.  If the plot of
#  crossvalidated sums of squares vs steps doesn not show a minimum, you
#  might have to increase num.steps.

devAskNewPage(ask = FALSE) 
hold = crossvalNparReg(formLatticeOutput,Z,
         PointPattern=grid2,M=0.5,num.steps = 200)
devAskNewPage(ask = TRUE)
         
#  Finally, createNparReg forms kernels using the lattice-based approach
#  and then performs kernel non-parametric regression.  Note that the
#  optimal crossvalidated smoothing parameter is in hold$k.  After that,
#  we plot the regression surface in a countour plot.

NparRegOut = createNparReg(formLatticeOutput,Z,PointPattern=grid2,k=hold$k)
plot(NparRegOut)
devAskNewPage(ask = FALSE)


