createNparReg <-
function(formLatticeOutput,Z,PointPattern=NULL,M=0.5,k,...)
{
#
#  This function only uses sparse matrices
#  It takes as input the formLatticeOutput and a column of responses
#  called Z.  Finally, k is the smoothing parameter and is the
#  number of steps in the random walk.
#
  if(class(formLatticeOutput)!="formLatticeOutput"){
       stop("Should be the output from the function formLattice")}
  if((M==0)|(M==1)){warning("Setting M to zero or one is ill-advised")}
  if(!is.null(PointPattern)){PointPattern <- as.matrix(PointPattern)}
  addQuantVarOutput = addQuantVar(formLatticeOutput = formLatticeOutput, Z=Z, 
  locations = PointPattern, will.plot=TRUE)
  p0 <- addQuantVarOutput$init.prob
  Z0 <- addQuantVarOutput$init.quantvar
  which.nodes <- addQuantVarOutput$which.nodes
  n <- length(PointPattern[,1])
  T<-makeTmatrix(formLatticeOutput=formLatticeOutput, M=M, sparse=TRUE)
  EW.locs <- formLatticeOutput$EW.locs
  NS.locs <- formLatticeOutput$NS.locs
  nodes <- formLatticeOutput$nodes
  pk <- Tkp(T,k=k,p=p0)
  zk <- Tkp(T,k=k,p=Z0)
  N <- length(NS.locs)*length(EW.locs)
  long <- rep(NA,N)
  zk[is.nan(zk)] <- mean(Z)
  pk[is.nan(pk)] <- 1
  long[as.numeric(rownames(nodes))] <- zk/pk
  NparRegOut = list(EW.locs = formLatticeOutput$EW.locs,
                    NS.locs = formLatticeOutput$NS.locs,
                    nodes = formLatticeOutput$nodes,
                    boundaryPoly = formLatticeOutput$poly,
                    hole.list = formLatticeOutput$hole.list,
                    PointPattern = PointPattern,
                    which.nodes = which.nodes,
                    NparRegDenom = zk,
                    NparRegNum = pk,
                    NparRegMap = long)
  class(NparRegOut) = "NparRegOut"
  return(NparRegOut)
}


