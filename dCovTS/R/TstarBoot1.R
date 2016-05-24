TstarBoot1 <- function(x,type,p,b,parallel=FALSE){  
 if(is.vector(x))stop('Multivariate time series only')
 if(!all(is.finite(x))) stop('Missing or infitive values')
 if (!is.numeric(x)) stop("'x' must be numeric") 
 n <- as.integer(NROW(x))
 q <- as.integer(NCOL(x))
 MaxLag <- n-2
 test <- function(j){
  kern <- kernelFun(type,j/p)
  if (kern==0){
   d=rep(0,b)
  } else {
  A <- crossDist(x,j)$A
  B <- crossDist(x,j)$B
  Atilde <- lapply(1:q, FUN=function(i) ATilde(A[[i]]))
  Btilde <- lapply(1:q, FUN=function(i) ATilde(B[[i]]))
  boot = function(Atilde,Btilde,j){
   Wtstar <- rbind(rnorm(n-j))
   Vrm <- matrix(NA,q,q)
    for (l in 1:q){
     for (f in 1:q){
       Vrm[l,f] <- sqrt((Wtstar%*%(Atilde[[l]]*Btilde[[f]])%*%t(Wtstar))/((n-j)^2))
     }
   }
   (n-j)*kern^2*sum(Vrm^2)
  }
  d=replicate(b,boot(Atilde,Btilde,j))
 }
 d
}
if(parallel==TRUE){
  closeAllConnections()
  cl <- makeCluster(detectCores())
  registerDoParallel(cl)
  clusterSetRNGStream(cl = cl, iseed = 9182)
  i <- seq_len(MaxLag)
  fe_call <- as.call(c(list (as.name("foreach"), i = i,.combine="+",.export=c("kernelFun","crossDist","ATilde")) ))
  fe <- eval(fe_call)
  res <- fe %dopar% test(i)
  stopCluster(cl)
}
else {
 res <- rowSums(sapply(1:MaxLag,function(i) test(i)))
}
return(res)
}


