RbootCV <- function(n,MaxLag,b,parallel=FALSE){ 
 x <- rnorm(n)
 if (missing(MaxLag) || MaxLag < 0) 
      stop("'MaxLag' must be greater than 1")
 X0 <- x[1:n]
 A0 <- as.matrix(dist(X0,diag=TRUE,upper=TRUE))
 total.meanA0 <- sum(A0)/(n^2)
 m1A0 <- sapply(1:NROW(A0), FUN=function(j) mean(A0[j,]))
 m2A0 <- sapply(1:NROW(A0), FUN=function(j) mean(A0[,j]))
 Atilde0 <- sapply(1:NROW(A0),1:NROW(A0),FUN=function(i,j) (A0[i,j]-m1A0[i]-m2A0[j]+total.meanA0)) 
 rstar <- function(k){
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
  Rstark = function(Atilde,Btilde,k){
   Wtstar <- rbind(rnorm(n-k))
   dcov <- sqrt((Wtstar%*%(Atilde*Btilde)%*%t(Wtstar))/((n-k)^2))
   dvarx <- sqrt(mean((Atilde0*Atilde0))*mean((Atilde0*Atilde0)))
   return(dcov/sqrt(dvarx))
  }
 return(replicate(b,Rstark(Atilde,Btilde,k)))
}
if(parallel==TRUE){
  closeAllConnections()
#  cl <- makeCluster(detectCores())
  cl <- makeCluster(2)
  registerDoParallel(cl)
  clusterSetRNGStream(cl = cl, iseed = 9182)
  i <- 1:MaxLag
  fe_call <- as.call( c(list (as.name("foreach"), i = i) ))
  fe <- eval(fe_call)
  Rstar <- fe %dopar% rstar(i)
  stopCluster(cl)
  quant <- sapply(1:MaxLag,FUN=function(j) quantile(Rstar[[j]],0.95))
  pv <- sapply(1:MaxLag,FUN=function(j) mean(Rstar[[j]]>=quant[j]))
  pvadj <- sapply(1:MaxLag,FUN=function(j) p.adjust(pv[j],method="fdr"))
  cvadj <- sapply(1:MaxLag,FUN=function(j) quantile(Rstar[[j]],1-pvadj[j]))
  res <- max(cvadj)
}
else {
 Rstar <- lapply(1:MaxLag,FUN=function(j) rstar(j))
 quant <- sapply(1:MaxLag,FUN=function(j) quantile(Rstar[[j]],0.95))
 pv <- sapply(1:MaxLag,FUN=function(j) mean(Rstar[[j]]>=quant[j+1]))
 pvadj <- sapply(1:MaxLag,FUN=function(j) p.adjust(pv[j],method="fdr"))
 cvadj <- sapply(1:MaxLag,FUN=function(j) quantile(Rstar[[j]],1-pvadj[j]))
 res <- max(cvadj)
}
return(res)
}




