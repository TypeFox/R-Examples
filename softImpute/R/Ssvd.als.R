Ssvd.als <-
  function (x, J = 2, thresh = 1e-05,lambda=0, maxit=100,trace.it=FALSE,warm.start=NULL,final.svd=TRUE) 
{
  this.call=match.call()
#  Compute the SVD via regularized ALS of a sparse matrix or SparsplusLowRank
  n <- dim(x)
  m <- n[2]
  n <- n[1]
  nz=n*m
  normx=norm(x,"F")# class the right method

   warm=FALSE
  if(!is.null(warm.start)){
    #must have u,D and V components
    if(!all(match(c("u","d","v"),names(warm.start),0)>0))stop("warm.start does not have components u, d and v")
    warm=TRUE
    D=warm.start$d
    JD=length(D)
    if(JD >= J){
      U=warm.start$u[,seq(J),drop=FALSE]
      V=warm.start$v[,seq(J),drop=FALSE]
      Dsq=D[seq(J)]
    }
    else{
      Dsq=c(D,rep(D[JD],J-JD))
      Ja=J-JD
      U=warm.start$u
      Ua=matrix(rnorm(n*Ja),n,Ja)
      Ua=Ua-U%*% (t(U)%*%Ua)
      Ua=svd(Ua)$u
      U=cbind(U,Ua)
      V=cbind(warm.start$v,matrix(0,m,Ja))
    }
  }
  else
    {
      V=matrix(0,m,J)
      U=matrix(rnorm(n*J),n,J)
      U=svd(U)$u
      Dsq=rep(1,J)# we call it Dsq because A=UD and B=VD and AB'=U Dsq V^T
    }

  ratio <- 1
  iter <- 0
  while ((ratio > thresh)&(iter<maxit)) {
    iter <- iter + 1
    U.old=U
    V.old=V
    Dsq.old=Dsq
   ## U step
    B=t(t(U)%*%x)
    if(lambda>0)B=UD(B,Dsq/(Dsq+lambda),m)
    Bsvd=svd(B)
    V=Bsvd$u
    Dsq=Bsvd$d
   # U=U%*%Bsvd$v
    ## V step
    obj=(.5*Frobsmlr(x,U,B,nx=normx)^2+lambda*sum(Dsq))/nz
    A=x%*%V
    if(lambda>0)A= UD(A,Dsq/(Dsq+lambda),n)
    Asvd=svd(A)
    U=Asvd$u
    Dsq=Asvd$d
    V=V %*% Asvd$v # just for computing the convergence criterion
    ratio=Frob(U.old,Dsq.old,V.old,U,Dsq,V)
    if(trace.it) cat("Ssvd.als:",iter, ":", "obj",format(round(obj,5)),"ratio", ratio, "\n")
  }
  if(iter==maxit)warning(paste("Convergence not achieved by",maxit,"iterations"))
  if((lambda>0)&&final.svd){
    ###Final cleanup of svd
    A=x%*%V
    Asvd=svd(A)
    U=Asvd$u
    Dsq=Asvd$d
    V=V %*% Asvd$v
    Dsq=pmax(Dsq-lambda,0)
  }
 ###

  out=list(u=U,d=Dsq,v=V)
  attr(out,"call")=this.call
  attr(out,"lambda")=lambda
  out
}
