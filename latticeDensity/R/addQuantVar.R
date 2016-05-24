addQuantVar <-
function(formLatticeOutput, Z, locations, will.plot=TRUE){
  #  
  #
  require(spdep)
  require(spatstat)
  if(class(formLatticeOutput)!="formLatticeOutput"){
       stop("Should be the output from the function formLattice")}
  nodes <- formLatticeOutput$nodes
  polyg <- formLatticeOutput$poly
  latt <- formLatticeOutput$latt
  n.observ <- nrow(locations)
  if(length(Z)!=n.observ){stop("The number of rows in the argument
  locations should be the same as the length of Z")}
  n.nodes <- nrow(nodes)
  Z <-  as.vector(Z)
#
#  Here we have to create ppp objects to use spatstat functions
#  to find nearest latts in nodes from locations
#
  temp <- bbox(rbind(locations,nodes))
  bound.vect <- c(temp[1,1],temp[1,2],temp[2,1],temp[2,2])
  X <- as.ppp(locations,W = bound.vect)
  Y <- as.ppp(nodes,W=bound.vect)
  closest <- nncross(X,Y)$which
#
#  Now I will plot the original locations
#  along with the nearest place in the nodes
#
  relocated.obj <- nodes[closest,]
  if(will.plot){
  plot(rbind(polyg,polyg[1,]),type="l")
  points(nodes)
  points(locations,pch=19,col=2)
  points(relocated.obj,pch=19,col=3,cex=0.5)
  }
#
#  The output will be a vector that gives an initial prob
#  at each node, depending on number of corresponding
#  locations
#
  temp = addObservations(formLatticeOutput,observations = locations)
  sums = function(x){sum(Z[closest==x])}
  out <- list(init.quantvar = as.numeric(lapply(1:n.nodes,FUN=sums))/n.observ,
             init.prob = temp$init.prob,
             which.nodes = closest)
  class(out) <- "addQuantVarOutput"
  return(out)
}

