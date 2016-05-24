bigtps <-
  function(x,y,nknots=NULL,nvec=NULL,rparm=NA,
           alpha=1,lambdas=NULL,se.fit=FALSE,
           rseed=1234,knotcheck=TRUE){
    ###### Fits Cubic Thin-Plate Spline
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: March 10, 2015
    
    ### initial info
    if(is.null(rseed)==FALSE){set.seed(rseed)}
    x <- as.matrix(x+0.0)
    n <- nrow(x)
    nx <- ncol(x)
    y <- as.matrix(y+0.0)
    yty <- sum(y^2)
    ysm <- sum(y)
    if(nrow(y)!=n){stop("Lengths of 'x' and 'y' must match.")}
    if(ncol(y)>1){stop("Response must be unidimensional (vector).")}
    if(nx>3){stop("Too many predictors. Use another function (bigssa or bigssp).")}
    nvec <- nvec[1]
    if(!is.null(nvec)){
      if(nvec<1L){
        if(nvec<=0){stop("Input 'nvec' is too small.")}
      } else{
        nvec <- as.integer(nvec)
        while(nvec<(nx+2)){nvec <- nvec+1L}
      }
    }
    
    ### check knots
    if(is.null(nknots)){
      if(nx==1L){nknots <- 30L} else if(nx==2L){nknots <- 100L} else {nknots <- 200L}
    }
    if(!is.null(rseed) && length(nknots)==1L){set.seed(rseed)}
    theknots <- NULL
    if(length(nknots)>1){
      kidx <- as.integer(nknots)
      nknots <- length(kidx)
      if(any(kidx<1L) || any(kidx>n)){stop("Input 'nknots' out of data index range.")}
    } else {
      kidx <- NA
    }
    
    ### round data
    if(!is.na(rparm[1])){
      xrng <- apply(x,2,range)
      xorig <- x
      yorig <- y
      gvec <- matrix(1L,n,1)
      kconst <- 1
      if(length(rparm) != nx){
        rparm <- rep(rparm[1],nx)
        rplog <- log(c(rparm[1],rparm[1]/2,rparm[1]/5),base=10)
        rpchk <- rep(FALSE,3)
        for(jj in 1:3){rpchk[jj] <- (rplog[jj]==as.integer(rplog[jj]))}
        if(!any(rpchk)){stop("Must set input 'rparm' such that rparm=a*(10^-b) with a in {1,2,5} and b>=1 (integer).")}
      } else{
        for(kk in 1:nx){
          rplog <- log(c(rparm[kk],rparm[kk]/2,rparm[kk]/5),base=10)
          rpchk <- rep(FALSE,3)
          for(jj in 1:3){rpchk[jj] <- (rplog[jj]==as.integer(rplog[jj]))}
          if(!any(rpchk)){stop("Must set input 'rparm' such that rparm=a*(10^-b) with a in {1,2,5} and b>=1 (integer).")}
        }
      }
      rxrng <- xrng
      for(j in 1:nx){
        if(rparm[j]<=0 || rparm[j]>=xrng[2,j]){stop("Must set input 'rparm' such that 0<rparm<max(x)")}
        #gvec = gvec + kconst*round((x[,j]-xrng[1,j])/rparm[j])
        #kconst = kconst*round(1+(xrng[2,j]-xrng[1,j])/rparm[j])
        rxrng[,j] <- round(rxrng[,j]/rparm[j])*rparm[j]
        gvec  <-  gvec + kconst*round((x[,j]-rxrng[1,j])/rparm[j])
        kconst  <-  kconst*round(1+(rxrng[2,j]-rxrng[1,j])/rparm[j])
        x[,j] <- as.matrix(round(x[,j]/rparm[j]))*rparm[j]
      }
      gvec <- as.factor(gvec)
      glindx <- split(cbind(1:n,y),gvec)
      if(!is.na(kidx[1])){theknots <- as.matrix(x[kidx,])}
      fs <- matrix(unlist(lapply(glindx,unifqsum)),ncol=3,byrow=TRUE)
      x <- as.matrix(x[fs[,1],])
      y <- fs[,2]
      w <- fs[,3]
      nunewr <- nrow(x)
    } else {
      nunewr <- n
      xorig <- yorig <- xrng <- NA
      w <- rep.int(1L,n)
      if(length(kidx)>1){theknots <- as.matrix(x[kidx,])}
    }
    
    ### get knot indices
    if(length(kidx)==1L){
      nknots <- as.integer(min(c(nunewr,nknots)))
      if(nknots<1){stop("Input 'nknots' must be positive integer.")}
      kidx <- sample(nunewr,nknots,prob=(w/n))
      theknots <- as.matrix(x[kidx,])
    }
    if(knotcheck){
      theknots <- unique(theknots)
      nknots <- nrow(theknots)
    }
    if(nknots<nx+2){stop("Input 'nknots' is too small. Need nknots-2 >= ncol(x).")}

    ### create penalty matrix
    if(nx<3){gconst <- 1} else{gconst <- (-1)}
    mqed <- gconst*(.Fortran("tpskersym",theknots,nknots,nx,matrix(0,nknots,nknots)))[[4]]
    sqed <- gconst*(.Fortran("tpsker",x,theknots,nunewr,nx,nknots,matrix(0,nunewr,nknots)))[[6]]
    
    ### check for eigenvalue decomposition of penalty matrix 
    if(is.null(nvec)){
      # tranforms problem to unconstrained
      Cqrd <- qr(cbind(1,theknots),LAPACK=TRUE)
      CqrQ <- qr.Q(Cqrd,complete=TRUE)
      CZ <- CqrQ[,(nx+2):nknots]
      rm(Cqrd,CqrQ)
      Qmat <- crossprod(CZ,mqed%*%CZ)
      # make design matrix and crossproduct matrix
      wsqrt <- sqrt(w)
      Kmat <- cbind(1,x)
      Jmat <- sqed%*%CZ
      rm(sqed)
      KtK <- crossprod(Kmat*wsqrt)
      KtJ <- crossprod(Kmat*w,Jmat)
      JtJ <- crossprod(Jmat*wsqrt)
      Kty <- crossprod(Kmat,y)
      Jty <- crossprod(Jmat,y)
      
    } else {
      # EVD of penalty matrix
      if(nvec<1L){
        qeig <- eigen(mqed,symmetric=TRUE)
        absval <- abs(qeig$val)
        qord <- order(absval,decreasing=TRUE)
        qvafs <- cumsum(absval[qord])/sum(absval)
        nvec <- which(qvafs>=nvec)[1]
        while(nvec<(nx+2)){nvec <- nvec+1L}
        Uvecs <- (qeig$vec[,qord])[,1:nvec]
        Uvals <- (qeig$val[qord])[1:nvec]
      } else{
        nvec <- as.integer(min(c(nvec[1],nknots)))
        if(nvec<1L){stop("Input 'nvec' must be positive integer.")}
        qeig <- eigen(mqed,symmetric=TRUE)
        qord <- order(abs(qeig$val),decreasing=TRUE)
        Uvecs <- (qeig$vec[,qord])[,1:nvec]
        Uvals <- (qeig$val[qord])[1:nvec]
      }
      rm(qeig,mqed)
      
      # tranforms problem to unconstrained
      Cqrd <- qr(crossprod(Uvecs,cbind(1,theknots)),LAPACK=TRUE)
      CqrQ <- qr.Q(Cqrd,complete=TRUE)
      CZ <- CqrQ[,(nx+2):nvec]
      rm(Cqrd,CqrQ)
      Qmat <- crossprod(CZ,diag(Uvals)%*%CZ)
      
      # make design matrix and crossproduct matrix
      wsqrt <- sqrt(w)
      Kmat <- cbind(1,x)
      Jmat <- sqed%*%(Uvecs%*%CZ)
      rm(sqed)
      KtK <- crossprod(Kmat*wsqrt)
      KtJ <- crossprod(Kmat*w,Jmat)
      JtJ <- crossprod(Jmat*wsqrt)
      Kty <- crossprod(Kmat,y)
      Jty <- crossprod(Jmat,y)
      
    } # end if(is.null(nvec))
    
    ### find optimal smoothing parameter 
    nbf <- nx+1L
    ncdim <- nrow(Qmat)
    if(is.null(lambdas)){
      nqmat <- matrix(0,nbf+ncdim,nbf+ncdim)
      nqmat[(nbf+1):(ncdim+nbf),(nbf+1):(ncdim+nbf)] <- n*Qmat
      newlam <- lamloop(10^-c(9:1),1,Kty,Jty,KtK,KtJ,JtJ,
                        Qmat,ncdim,n,alpha,yty,nbf)
      gcvopt <- nlm(f=gcvcss,p=log(newlam),yty=yty,
                    xtx=rbind(cbind(KtK,KtJ),cbind(t(KtJ),JtJ)),
                    xty=rbind(Kty,Jty),nqmat=nqmat,ndpts=n,alpha=alpha)
      lambda <- exp(gcvopt$est)
    } else if(length(lambdas)>1){
      if(any(lambdas<0)){stop("Input 'lambdas' must be nonnegative.")}
      lambda <- lamloop(lambdas,1,Kty,Jty,KtK,KtJ,JtJ,
                        Qmat,ncdim,n,alpha,yty,nbf)
    } else {
      lambda <- lambdas[1]
      if(lambda<0){stop("Input 'lambdas' must be nonnegative.")}
    }
    
    ### get final estimates
    fxhat <- lamcoef(lambda,1,Kty,Jty,KtK,KtJ,JtJ,Qmat,ncdim,n,alpha,yty,nbf)
    fhat <- fxhat[[1]]
    dchat <- fhat[1:(nbf+ncdim)]
    yhat <- cbind(Kmat,Jmat)%*%dchat
    sseval <- fhat[nbf+ncdim+1]
    effdf <- fhat[nbf+ncdim+2]
    mevar <- sseval/(n-effdf)
    gcv <- n*sseval/((n-alpha*effdf)^2)
    aic <- n*(1+log(2*pi)) + n*log(sseval/n) + effdf*2
    bic <- n*(1+log(2*pi)) + n*log(sseval/n) + effdf*log(n)
    csqrt <- sqrt(mevar)*fxhat[[2]]
    
    ### posterior variance
    pse <- NA
    if(se.fit){pse <- sqrt(postvar(Kmat,Jmat,csqrt))}
    if(is.null(nvec)){
      dchat <- c(dchat[1:(nx+1)],CZ%*%dchat[(nx+2):(ncdim+nbf)])
      csqrt <- rbind(csqrt[1:(nx+1),],CZ%*%csqrt[(nx+2):(ncdim+nbf),])
    } else {
      dchat <- c(dchat[1:(nx+1)],Uvecs%*%CZ%*%dchat[(nx+2):(ncdim+nbf)])
      csqrt <- rbind(csqrt[1:(nx+1),],Uvecs%*%CZ%*%csqrt[(nx+2):(ncdim+nbf),])
    }
    
    ### calculate vaf
    mval <- ysm/n
    vaf <- 1 - sseval/(yty-n*(mval^2))
    
    ### collect results
    if(is.na(rparm[1])==FALSE){
      xunique <- x
      yunique <- y/w
      y <- yorig
      x <- xorig 
    } else {xunique <- yunique <- w <- NA}
    ndf <- data.frame(n=n,df=effdf,row.names="")
    tpsfit <- list(fitted.values=yhat,se.fit=pse,x=x,y=y,
                   xunique=xunique,yunique=yunique,funique=w,sigma=sqrt(mevar),
                   ndf=ndf,info=c(gcv=gcv,rsq=vaf,aic=aic,bic=bic),
                   myknots=theknots,nvec=nvec,rparm=rparm,
                   lambda=lambda,coef=dchat,coef.csqrt=csqrt)
    class(tpsfit) <- "bigtps"
    return(tpsfit)
    
  }