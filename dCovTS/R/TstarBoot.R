TstarBoot <- function(x,type,testType,p,b,parallel=FALSE){ 
 n <- length(x)
 if (is.matrix(x)) {
  if (!NCOL(x)==1) stop('Univariate time series only')
 } else {
  x <- c(x)
 }
 if(!all(is.finite(x))) stop('Missing or infitive values')
 if (!is.numeric(x)) 
     stop("'x' must be numeric") 
 if((b==0) | missing(b)) stop('b must be grater than 0')
 X0 <- x[1:n]
 A0 <- as.matrix(dist(X0,diag=TRUE,upper=TRUE))
 total.meanA0 <- sum(A0)/(n^2)
 m1A0 <- sapply(1:NROW(A0), FUN=function(j) mean(A0[j,]))
 m2A0 <- sapply(1:NROW(A0), FUN=function(j) mean(A0[,j]))
 Atilde0 <- sapply(1:NROW(A0),1:NROW(A0),FUN=function(i,j) (A0[i,j]-m1A0[i]-m2A0[j]+total.meanA0)) 
 MaxLag <- n-1
 test <- function(k){
  kern <- kernelFun(type,k/p)
  if (kern==0){
   d=rep(0,b)
  } else {
  Y <- x[1:(n-k)]
  X <- x[(1+k):n]
  A <- as.matrix(dist(X,diag=TRUE,upper=TRUE))
  B <- as.matrix(dist(Y,diag=TRUE,upper=TRUE))
  total.meanA <- sum(A)/(NROW(A)^2)
  m1A <- sapply(1:NROW(A), FUN=function(j) mean(A[j,]))
  m2A <- sapply(1:NROW(A), FUN=function(j) mean(A[,j]))
  Atilde <- sapply(1:NROW(A),1:NROW(A),FUN=function(i,j) (A[i,j]-m1A[i]-m2A[j]+total.meanA))
  total.meanB <- sum(B)/(NROW(B)^2)
  m1B <- sapply(1:NROW(B), FUN=function(j) mean(B[j,]))
  m2B <- sapply(1:NROW(B), FUN=function(j) mean(B[,j]))
  Btilde <- sapply(1:NROW(B),1:NROW(B),FUN=function(i,j) (B[i,j]-m1B[i]-m2B[j]+total.meanB))  
  bootCov = function(Atilde,Btilde,k){
   Wtstar <- rbind(rnorm(n-k))
   V <- sqrt((Wtstar%*%(Atilde*Btilde)%*%t(Wtstar))/((n-k)^2))
   return((n-k)*kern^2*V^2)
  }
  bootCor = function(Atilde,Btilde,k){
    Wtstar <- rbind(rnorm(n-k))
    dcov <- sqrt((Wtstar%*%(Atilde*Btilde)%*%t(Wtstar))/((n-k)^2))
    dvarx <- sqrt(mean((Atilde0*Atilde0))*mean((Atilde0*Atilde0)))
    V <- dcov/sqrt(dvarx)
    return((n-k)*kern^2*V^2)
  } 
  if (testType=="covariance"){
   d=replicate(b,bootCov(Atilde,Btilde,k))
  } else {
   d=replicate(b,bootCor(Atilde,Btilde,k))
  }
 }
d
}
if(parallel==TRUE){
  closeAllConnections()
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  clusterSetRNGStream(cl = cl, iseed = 9182)
  i <- seq_len(MaxLag)
  fe_call <- as.call( c(list (as.name("foreach"), i = i,.combine="+",.export="kernelFun") ))
  fe <- eval(fe_call)
  res <- fe %dopar% test(i)
  stopCluster(cl)
}
else {
 res <- rowSums(sapply(1:MaxLag,function(k) test(k)))
}
return(res)
}
