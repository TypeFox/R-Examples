createDensity <-
function(formLatticeOutput,PointPattern=NULL,M=0.5,k,intensity=FALSE, 
sparse=TRUE, ...)
{
# 
  if(sparse){require(spam)}
  if(class(formLatticeOutput)!="formLatticeOutput"){
       stop("Should be the output from the function formLattice")}
  if((M==0)|(M==1)){warning("Setting M to zero or one is ill-advised")}
  init.prob = addObservations(formLatticeOutput, PointPattern)
  p0 <- as.vector(init.prob$init.prob)
  if(!is.null(PointPattern)){PointPattern <- as.matrix(PointPattern)}
  poly.area <- areaRegion(formLatticeOutput)
  n <- length(PointPattern[,1])
  T<-makeTmatrix(formLatticeOutput, M=M, sparse=sparse)
  EW.locs <- formLatticeOutput$EW.locs
  NS.locs <- formLatticeOutput$NS.locs
  nodes <- formLatticeOutput$nodes
  z <- Tkp(T,k=k,p=p0)
  N <- length(NS.locs)*length(EW.locs)
  long <- rep(NA,N)
  if(intensity){
  long[as.numeric(rownames(nodes))] <- 
       n*z*length(z)/poly.area
  }else{
  long[as.numeric(rownames(nodes))] <- 
       z*length(z)/poly.area}
  densityOut     <- list(EW.locs = formLatticeOutput$EW.locs,
                    NS.locs = formLatticeOutput$NS.locs,
                    nodes = formLatticeOutput$nodes,
                    boundaryPoly = formLatticeOutput$poly,
                    hole.list = formLatticeOutput$hole.list,
                    PointPattern = PointPattern,
                    probs = z,
                    densityOut = long,
                    area = poly.area)
  class(densityOut) <- "densityOut"
  return(densityOut)
}


