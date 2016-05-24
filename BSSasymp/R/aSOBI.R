aSOBI <- function(X, k=12, a=4, eps = 1e-06, maxiter = 1000)
{
  if (length(k)==1) k <- 1:k  
  n <- dim(X)[1]
  p <- dim(X)[2]   
  K <- length(k)
     
  sobi <- SOBI(X,k=k,eps=eps,maxiter=maxiter) 
  W0 <- sobi$W
  Z <- sobi$S

  CM <- array(0,dim=c(p,p,K))
  for(i in 1:K){
    S <- crossprod(Z[1:(n-k[i]),],Z[(k[i]+1):n,])/(n-k[i])
    CM[,,i] <- (1/2)*(S+t(S))
  } 
  
  V <- diag(p) 
  V0 <- V
  W <- matrix(0,p,p)
  iter <- 0
 
  while (TRUE){
    iter <- iter+1
    W <- matrix(0,p,p)
    for(mi in 1:K){
      W <- W+2*tcrossprod(V,CM[,,mi])*matrix(rep(diag(g(tcrossprod(tcrossprod(V,CM[,,mi]), V),a)),p),p,p)
    }  
    V <- crossprod(solve(mat.sqrt(tcrossprod(W,W))),W)
    if(mat.norm(V-V0)<eps) break
    if(iter==maxiter) stop("maxiter reached without convergence") 
    V0<-V
  }
  W <- V%*%W0
  S <- tcrossprod(Z,V)
    
  acs <- acf(S, lag.max=max(k), plot=FALSE)$acf
  ssq_ac <- NULL
  for(j in 1:p){
    ssq_ac[j]<-sum(abs(acs[k+1,j,j])^a)
  }
  ord <- order(ssq_ac, decreasing=TRUE)
  P <- matrix(0,p,p)
  for(j in 1:p){
    P[j,ord[j]]<-1
  }  
  W <- P%*%W
  W <- diag(sign(rowMeans(W)))%*%W

  S <- tcrossprod(sweep(X,2,colMeans(X),"-"),W)
  S <- ts(S, names=paste("Series",1:p))

  res <- list(W=W, S=S, k=k, a=a) 
  class(res) <- "bss"
  res
}
 
