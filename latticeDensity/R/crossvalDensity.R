crossvalDensity <- 
function(formLatticeOutput,PointPattern, M=0.5, num.steps = 200, sparse=TRUE){
#
#
#
  if(class(formLatticeOutput)!="formLatticeOutput"){
       stop("Should be the output from the function formLattice")}
  if((M==0)|(M==1)){warning("Setting M to zero or one is ill-advised")}
  PointPattern <- as.matrix(PointPattern)
  addObsOutput <- addObservations(formLatticeOutput,PointPattern)
  p0 <- addObsOutput$init.prob
  which.nodes <- addObsOutput$which.nodes
  NN <- length(p0)
  n <- length(PointPattern[,1])
  poly.area <- areaRegion(formLatticeOutput)
#
if(sparse){
  #
  hold.del.prob <- matrix(nrow=NN, ncol=n, NA)
  for (i in 1:n){
   temp <- addObservations(formLatticeOutput, PointPattern[-i,], will.plot=FALSE)
   hold.del.prob[,i] <- temp$init.prob
  }
  T <- makeTmatrix(formLatticeOutput,M = M, sparse=TRUE)
  Tkinit <-  p0
  first.term <- rep(NA,num.steps)
  second.term <- rep(NA,num.steps)
  for(k in 1:num.steps){
    Tkinit <- T%*%Tkinit
    first.term[k] <- (NN/poly.area)*sum(Tkinit*Tkinit)
    hold.del.prob <- T%*%hold.del.prob
    temp = sum(hold.del.prob[cbind(which.nodes,1:n)])
    second.term[k] <- (NN/poly.area)*(2/n)*temp
} 
#
} else {
#
  T = makeTmatrix(formLatticeOutput,M = 0.5, sparse=FALSE)
  Tnew <- diag(x=1,nrow=NN,ncol=NN)
  first.term <- rep(NA,num.steps)
  second.term <- rep(NA,num.steps)
  for(k in 1:num.steps){
    Tnew <- Tnew%*%T
    probs <- Tnew%*%p0
    first.term[k] <- (NN/poly.area)*sum(probs^2)
    second.term[k] <- (NN/poly.area)*(2/n)*((n^2/(n-1))*p0%*%probs - (n/(n-1))*sum(diag(Tnew)*p0))
  } 
 }
  ucv <- first.term - second.term
  plot(ucv,type="l")
  k <- which.min(ucv)
  out <- list(ucv=ucv,k=k)
  return(out)
  }











