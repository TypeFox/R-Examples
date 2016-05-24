addObservations <-
function(formLatticeOutput, observations, will.plot=TRUE){
  #  Version 20July.1
  #  Changed output to only init.prob
  #
  require(spdep)
  require(spatstat)
  if(class(formLatticeOutput)!="formLatticeOutput"){
       stop("Should be the output from the function formLattice")}
  nodes <- formLatticeOutput$nodes
  polyg <- formLatticeOutput$poly
  latt <- formLatticeOutput$latt
  n.observ <- nrow(observations)
  n.nodes <- nrow(nodes)
#
#  Here we have to create ppp objects to use spatstat functions
#  to find nearest latts in nodes from observations
#
  temp <- bbox(rbind(observations,nodes))
  bound.vect <- c(temp[1,1],temp[1,2],temp[2,1],temp[2,2])
  X <- as.ppp(observations,W = bound.vect)
  Y <- as.ppp(nodes,W=bound.vect)
  closest <- nncross(X,Y)$which
#
#  Now I will plot the original observations
#  along with the nearest place in the nodes
#
  relocated.obj <- nodes[closest,]
  if(will.plot){
  plot(rbind(polyg,polyg[1,]),type="l")
  points(nodes)
  points(observations,pch=19,col=2)
  points(relocated.obj,pch=19,col=3,cex=0.5)
  }
#
#  The output will be a vector that gives an initial prob
#  at each node, depending on number of corresponding
#  observations
#
  out <- list(init.prob = tabulate(closest,nbins=n.nodes)/n.observ,
             which.nodes = closest)
  class(out) <- "initProbObject"
  return(out)
}

