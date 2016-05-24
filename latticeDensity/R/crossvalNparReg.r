crossvalNparReg <- 
function(formLatticeOutput,Z,PointPattern,M=0.5,num.steps = 200,plots=TRUE){
#
#                                                             
#
  if(class(formLatticeOutput)!="formLatticeOutput"){
       stop("Should be the output from the function formLattice")}
  if((M==0)|(M==1)){warning("Setting M to zero or one is ill-advised")}
  PointPattern <- as.matrix(PointPattern)
  addQuantVarOut <- addQuantVar(formLatticeOutput, Z=Z, 
      locations = PointPattern, will.plot=plots)
  p0 <- addQuantVarOut$init.prob
# Z0 <- addQuantVarOut$init.quantvar
  which.nodes <- addQuantVarOut$which.nodes
  NN <- length(p0)
  n <- length(PointPattern[,1])
  hold.del.prob <- matrix(nrow=NN, ncol=n, NA)
  hold.del.Z <- matrix(nrow=NN, ncol=n, NA)
  for (i in 1:n){
    temp <- addQuantVar(formLatticeOutput, Z=Z[-i], 
        locations = PointPattern[-i,], will.plot=plots)
    hold.del.prob[,i] <- temp$init.prob
    hold.del.Z[,i] <- temp$init.quantvar
  }
  T <- makeTmatrix(formLatticeOutput,M = M, sparse=TRUE)
# # Tkinit <-  p0
# # TZinit <-  Z0
  SumSq <- rep(NA,num.steps)
  for(k in 1:num.steps){
  #  Tkinit <- T%*%Tkinit
  #  TZinit <- T%*%TZinit
    hold.del.prob <- T%*%hold.del.prob
    hold.del.Z <- T%*%hold.del.Z
    predictions <- diag((hold.del.Z/hold.del.prob)[which.nodes,])
    nnan <- sum(is.nan(predictions))
    predictions[is.nan(predictions)] <- (rep(sum(Z),nnan) - 
             Z[is.nan(predictions)])/(n-1)
    deleted.residuals <- predictions - Z
    SumSq[k] <- sum((deleted.residuals)^2)
   }
   plot(SumSq,type="l")
   k <- which.min(SumSq)
   out <- list(SumSq=SumSq,k=k)
   return(out) 
}











