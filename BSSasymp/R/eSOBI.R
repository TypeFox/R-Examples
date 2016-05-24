eSOBI <- function(X, taus=taus_def, M=200, fast=TRUE, eps=1e-06, maxiter=1000)
{
  n <- dim(X)[1]
  p <- dim(X)[2]   
  K <- length(taus)
  Ws <- array(0,c(p,p,K))
  sum_var <- rep(Inf,K)

  if(fast){
    Z<-SOBI(X,taus[[1]],eps=eps,maxiter=maxiter)$S
   
    F_tau<-array(0,c(p,p,M+1))
    for(m in 0:M){
      F_tau[,,m+1]<-tcrossprod(t(Z[1:(n-m),]),t(Z[(m+1):n,]))/(n-m)
    }
    Beta<-2*diag(p)+matrix(1,p,p)

    sasv<-NULL
    for(i in 1:length(taus)){
      Ki <- length(taus[[i]])     
      Lambda<-array(0,c(p,p,Ki))
      for(k in 1:Ki){
        for(j in 1:p){  
          Lambda[j,j,k]<-F_tau[j,j,taus[[i]][k]+1] 
        }
      }
      sum_var[i]<- sum(.C("ascov_all", as.double(as.vector(F_tau)),as.double(as.vector(Lambda)),as.double(taus[[i]]),as.integer(c(p,M,Ki)),as.double(as.vector(Beta)),res=double(p*(p-1)), PACKAGE="BSSasymp")$res)
    }
    b<-which.min(sum_var)
    W<-SOBI(X,taus[[b]],eps=eps,maxiter=maxiter)$W
  }else{ 
    for(i in 1:K){ 
      sobi <- tryCatch(SOBI(X, taus[[i]], eps=eps, maxiter=maxiter), error=function(e) 0)
      if(is.list(sobi)){
        Ws[,,i] <- sobi$W
        ascov <- ASCOV_SOBI_estN(sobi$S,taus[[i]], mixed=FALSE, M=M)$COV_W*n
        sum_var[i] <- sum(diag(ascov))
      }
    }
    b <- which.min(sum_var)
    W <- Ws[,,b]
  }  

  S <- tcrossprod(X,W)
  S <- sweep(S,2,colMeans(S),"-") 
  res <- list(W=W, S=S, taus_used=taus[[b]], sum_var=sum_var)
  class(res) <- "bss"
  res
}

