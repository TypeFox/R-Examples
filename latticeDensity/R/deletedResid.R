deletedResid <- 
function(formLatticeOutput,Z,PointPattern,M=0.5,k){
#
#                                                            
#
  if(class(formLatticeOutput)!="formLatticeOutput"){
       stop("Should be the output from the function formLattice")}
  if((M==0)|(M==1)){warning("Setting M to zero or one is ill-advised")}
  PointPattern <- as.matrix(PointPattern)
  addQuantVarOut <- addQuantVar(formLatticeOutput, Z=Z, 
      locations = PointPattern, will.plot=TRUE)
  pk <- addQuantVarOut$init.prob
  Zk <- addQuantVarOut$init.quantvar
  which.nodes <- addQuantVarOut$which.nodes
  NN <- length(pk)
  n <- length(PointPattern[,1])
  hold.del.prob <- matrix(nrow=NN, ncol=n, NA)
  hold.del.Z <- matrix(nrow=NN, ncol=n, NA)
  for (i in 1:n){
    temp <- addQuantVar(formLatticeOutput, Z=Z[-i], 
        locations = PointPattern[-i,], will.plot=TRUE)
    hold.del.prob[,i] <- temp$init.prob
    hold.del.Z[,i] <- temp$init.quantvar
  }
  T <- makeTmatrix(formLatticeOutput,M = M, sparse=TRUE)
  for(k in 1:k){
    pk <- T%*%pk
    Zk <- T%*%Zk
    hold.del.prob <- T%*%hold.del.prob
    hold.del.Z <- T%*%hold.del.Z
}
    fitted <- (Zk/pk)[which.nodes,]
    fitted[is.nan(fitted)] <- mean(Z)
    residuals <- fitted - Z
    del.predictions <- diag((hold.del.Z/hold.del.prob)[which.nodes,])
    nnan <- sum(is.nan(del.predictions))
    del.predictions[is.nan(del.predictions)] <- (rep(sum(Z),nnan) - 
             Z[is.nan(del.predictions)])/(n-1)
    deleted.residuals <- del.predictions - Z
   return(list(deleted.residuals=deleted.residuals,residuals=residuals)) 
}











