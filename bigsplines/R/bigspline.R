bigspline <-
  function(x,y,type="cub",nknots=30,rparm=0.01,xmin=min(x),
           xmax=max(x),alpha=1,lambdas=NULL,se.fit=FALSE,
           rseed=1234,knotcheck=TRUE){
    ###### Fits Smoothing Spline
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: January 14, 2016
    
    ### initial info
    if(!is.null(rseed)){set.seed(rseed)}
    x <- as.matrix(x+0.0)
    n <- nrow(x)
    nx <- ncol(x)
    y <- as.matrix(y+0.0)
    yty <- sum(y^2)
    ysm <- sum(y)
    if(nrow(y) != n){stop("Lengths of 'x' and 'y' must match.")}
    if(ncol(y) > 1){stop("Response must be unidimensional (vector).")}
    nknots <- as.integer(nknots[1])
    if(nknots < 1L){stop("Input 'nknots' must be positive integer.")}
    type <- type[1]
    if(!any(type==c("cub","cub0","per","lin"))){stop("Must set type to 'cub', 'cub0', 'per', or 'lin'.")}
    
    ### check inputs and transform data
    if(nx > 1){stop("Too many predictors. Use another function (bigssa or bigssp).")}
    if(xmin > xmax){stop("Input 'xmin' must be less than input 'xmax'.")}
    xrng <- matrix(c(xmin,xmax),2,1)
    x <- (x-xrng[1])/(xrng[2]-xrng[1])
    
    ### round data
    rparm <- rparm[1]
    if(!is.na(rparm)){
      if(rparm<=0 || rparm>=1){stop("Must set input 'rparm' such that 0<rparm<1")}
      rplog <- log(c(rparm,rparm/2,rparm/5),base=10)
      rpchk <- rep(FALSE,3)
      for(jj in 1:3){rpchk[jj] <- (rplog[jj]==as.integer(rplog[jj]))}
      if(!any(rpchk)){stop("Must set input 'rparm' such that rparm=a*(10^-b) with a in {1,2,5} and b>=1 (integer).")}
      xorig <- x
      yorig <- y
      gidx <- round(x/rparm)
      gvec <- as.integer(1L+gidx)
      x <- as.matrix(seq(0,1,by=rparm))
      nunewr <- nrow(x)
      yw <- (.Fortran("sumfreq",y,gvec,n,nunewr,double(nunewr),
                      integer(nunewr),PACKAGE="bigsplines"))[5:6]
      widx <- which(yw[[2]]>0)
      y <- as.matrix(yw[[1]][widx])
      w <- yw[[2]][widx]
      rm(yw)
      x <- as.matrix(x[widx])
      nunewr <- nrow(x)
    } else {
      nunewr <- n
      xorig <- yorig <- NA
      w <- 1
    }
    
    ### get knots
    kidx <- binsamp(x,matrix(c(0,1),2,1),min(nknots,nrow(x)),1L)
    theknots <- as.matrix(x[kidx])
    nknots <- length(kidx)
    if(knotcheck){
      theknots <- unique(theknots)
      nknots <- nrow(theknots)
    }
    
    ### make marginal reproducing kernel matrices
    if(type=="cub"){
      Kmat <- cbind(1,x-0.5)
      nbf <- 2L
      Jmat <- (.Fortran("cubker",x,theknots,nunewr,nknots,
                        matrix(0,nunewr,nknots),PACKAGE="bigsplines"))[[5]]
      Qmat <- (.Fortran("cubkersym",theknots,nknots,
                        matrix(0,nknots,nknots),PACKAGE="bigsplines"))[[3]]
    } else if(type=="per"){
      Kmat <- matrix(1,nunewr)
      nbf <- 1L
      Jmat <- (.Fortran("perker",x,theknots,nunewr,nknots,
                        matrix(0,nunewr,nknots),PACKAGE="bigsplines"))[[5]]
      Qmat <- (.Fortran("perkersym",theknots,nknots,
                        matrix(0,nknots,nknots),PACKAGE="bigsplines"))[[3]]
    } else if(type=="lin"){
      Kmat <- matrix(1,nunewr)
      nbf <- 1L
      Jmat <- (.Fortran("linker",x,theknots,nunewr,nknots,
                        matrix(0,nunewr,nknots),PACKAGE="bigsplines"))[[5]]
      Qmat <- (.Fortran("linkersym",theknots,nknots,
                        matrix(0,nknots,nknots),PACKAGE="bigsplines"))[[3]]
    } else {
      Kmat <- cbind(1,x)
      nbf <- 2L
      Jmat <- (.Fortran("cubkerz",x,theknots,nunewr,nknots,
                        matrix(0,nunewr,nknots),PACKAGE="bigsplines"))[[5]]
      Qmat <- (.Fortran("cubkerzsym",theknots,nknots,
                        matrix(0,nknots,nknots),PACKAGE="bigsplines"))[[3]]
    }
    
    ### get cross-product matrices
    wsqrt <- sqrt(w)
    KtJ <- crossprod(Kmat*w,Jmat)
    KtK <- crossprod(Kmat*wsqrt)
    JtJ <- crossprod(Jmat*wsqrt)
    Kty <- crossprod(Kmat,y)
    Jty <- crossprod(Jmat,y)
    
    ### find optimal smoothing parameter 
    if(is.null(lambdas)){
      nqmat <- matrix(0,nbf+nknots,nbf+nknots)
      nqmat[(nbf+1):(nknots+nbf),(nbf+1):(nknots+nbf)] <- n*Qmat
      newlam <- lamloop(10^-c(9:1),1,Kty,Jty,KtK,KtJ,JtJ,
                        Qmat,nknots,n,alpha,yty,nbf)
      gcvopt <- nlm(f=gcvcss,p=log(newlam),yty=yty,
                    xtx=rbind(cbind(KtK,KtJ),cbind(t(KtJ),JtJ)),
                    xty=rbind(Kty,Jty),nqmat=nqmat,ndpts=n,alpha=alpha)
      lambda <- exp(gcvopt$est)
    } else if(length(lambdas)>1){
      if(any(lambdas<0)){stop("Input 'lambdas' must be nonnegative.")}
      lambda <- lamloop(lambdas,1,Kty,Jty,KtK,KtJ,JtJ,
                        Qmat,nknots,n,alpha,yty,nbf)
    } else {
      lambda <- lambdas[1]
      if(lambda<0){stop("Input 'lambdas' must be nonnegative.")}
    }
    
    ### get final estimates
    fxhat <- lamcoef(lambda,1,Kty,Jty,KtK,KtJ,JtJ,
                     Qmat,nknots,n,alpha,yty,nbf)
    fhat <- fxhat[[1]]
    dchat <- fhat[1:(nbf+nknots)]
    yhat <- cbind(Kmat,Jmat)%*%dchat
    sseval <- fhat[nbf+nknots+1]
    effdf <- fhat[nbf+nknots+2]
    mevar <- sseval/(n-effdf)
    gcv <- n*sseval/((n-alpha*effdf)^2)
    aic <- n*(1+log(2*pi)) + n*log(sseval/n) + effdf*2
    bic <- n*(1+log(2*pi)) + n*log(sseval/n) + effdf*log(n)
    csqrt <- sqrt(mevar)*fxhat[[2]]
    
    ### posterior variance
    pse <- NA
    if(se.fit){pse <- sqrt(postvar(Kmat,Jmat,csqrt))}
    
    ### calculate vaf
    mval <- ysm/n
    vaf <- 1 - sseval/(yty-n*(mval^2))
    
    ### collect results
    x <- x*(xrng[2]-xrng[1]) + xrng[1]
    if(!is.na(rparm[1])){
      xunique <- x
      yunique <- y/w
      y <- yorig
      x <- xorig*(xrng[2]-xrng[1]) + xrng[1]
    } else {
      xunique <- yunique <- w <- NA
    }
    ndf <- data.frame(n=n,df=effdf,row.names="")
    sspfit <- list(fitted.values=yhat,se.fit=pse,x=x,y=y,type=type,
                   xunique=xunique,yunique=yunique,funique=w,sigma=sqrt(mevar),
                   ndf=ndf,info=c(gcv=gcv,rsq=vaf,aic=aic,bic=bic),
                   xrng=xrng,myknots=theknots,rparm=rparm,lambda=lambda,
                   coef=dchat,coef.csqrt=csqrt)
    class(sspfit) <- "bigspline"
    return(sspfit)
    
  }