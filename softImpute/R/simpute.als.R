simpute.als<-function (x, J = 2, thresh = 1e-05,lambda=0,maxit=100,trace.it=TRUE, warm.start=NULL, final.svd=TRUE) 
{### x is a matrix, possibly with NAs
    n <- dim(x)
    m <- n[2]
    n <- n[1]
    this.call=match.call()
    a=names(attributes(x))
    binames=c("biScale:row","biScale:column")
    if(all(match(binames,a,FALSE))){
      biats=attributes(x)[binames]
    } else biats=NULL

    xnas <- is.na(x)
    nz=m*n-sum(xnas)
    xfill <- x
      if(!is.null(warm.start)){
    #must have u,d and v components
    if(!all(match(c("u","d","v"),names(warm.start),0)>0))stop("warm.start does not have components u, d and v")
    warm=TRUE
    D=warm.start$d
    JD=sum(D>0)
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
  xfill[xnas]=(U%*%(Dsq*t(V)))[xnas]  
  }
  else
    {
      V=matrix(0,m,J)
      U=matrix(rnorm(n*J),n,J)
      U=svd(U)$u
      Dsq=rep(1,J)# we call it Dsq because A=UD and B=VD and AB'=U Dsq V^T
      xfill[xnas]=0
    }
    ratio <- 1
    iter <- 0
   while ((ratio > thresh)&(iter<maxit)) {
        iter <- iter + 1
        U.old=U
        V.old=V
        Dsq.old=Dsq
       ## U step
        B=t(U)%*%xfill
        if(lambda>0)B=B*(Dsq/(Dsq+lambda))
        Bsvd=svd(t(B))
        V=Bsvd$u
        Dsq=(Bsvd$d)
        U=U%*%Bsvd$v
        xhat=U %*%(Dsq*t(V))
        xfill[xnas]=xhat[xnas]
 ###The next line we could have done later; this is to match with sparse version
       if(trace.it) obj=(.5*sum( (xfill-xhat)[!xnas]^2)+lambda*sum(Dsq))/nz
        ## V step
        A=t(xfill%*%V)
        if(lambda>0)A=A*(Dsq/(Dsq+lambda))
        Asvd=svd(t(A))
        U=Asvd$u
        Dsq=Asvd$d
        V=V %*% Asvd$v
        xhat=U %*%(Dsq*t(V))
        xfill[xnas]=xhat[xnas]
        ratio=Frob(U.old,Dsq.old,V.old,U,Dsq,V)
        if(trace.it) cat(iter, ":", "obj",format(round(obj,5)),"ratio", ratio, "\n")
    }
    if(iter==maxit)warning(paste("Convergence not achieved by",maxit,"iterations"))

    if(lambda>0&final.svd){
      U=xfill%*%V
      sU=svd(U)
      U=sU$u
      Dsq=sU$d
      V=V%*%sU$v
      Dsq=pmax(Dsq-lambda,0)
      if(trace.it){
        xhat=U %*%(Dsq*t(V))
        obj=(.5*sum( (xfill-xhat)[!xnas]^2)+lambda*sum(Dsq))/nz
        cat("final SVD:", "obj",format(round(obj,5)),"\n")
    }

  }
 J=min(sum(Dsq>0)+1,J)
    out=list(u=U[,seq(J)],d=Dsq[seq(J)],v=V[,seq(J)])
  attributes(out)=c(attributes(out),list(lambda=lambda,call=this.call),biats)
  out

}

