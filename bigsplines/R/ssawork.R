ssawork <-
  function(formula,ssamk){
    
    ### check inputs
    if(formula==ssamk$formula){
      Etab <- ssamk$et
      mfdim <- dim(Etab)
      tnames <- colnames(Etab)
      oldnames <- rownames(ssamk$et)[2:(ssamk$nxvar+1L)]
    } else{
      
      # get formula
      trmfrm <- terms.formula(formula)
      ychk <- attr(trmfrm,"response")
      if(ychk!=1L){stop("Must include same response in formula when entering 'makessa' object.")}
      et <- attr(trmfrm,"factors")
      mfdim <- dim(et)
      newnames <- rownames(et)
      tnames <- colnames(et)
      mtidx <- match(tnames,colnames(ssamk$et))
      if(any(is.na(mtidx))){
        idxna <- which(is.na(mtidx))
        for(j in idxna){
          newsplit <- unlist(strsplit(tnames[j],":"))
          lnwsp <- length(newsplit)
          if(lnwsp==1L){
            stop("Cannot add new effects to formula when entering 'makessa' object.")
          } else if(lnwsp==2L){
            newmatch <- match(paste(newsplit[2],newsplit[1],sep=":"),colnames(ssamk$et))
            if(is.na(newmatch)){
              stop("Cannot add new effects to formula when entering 'makessa' object.")
            } else{
              tnames[j] <- paste(newsplit[2],newsplit[1],sep=":")
            }
          } else if(lnwsp==3L){
            myperms <- matrix(c(3,2,3,1,2,2,3,1,3,1,1,1,2,2,3),ncol=3)
            newmatch <- rep(NA,5)
            for(k in 1:5){newmatch[k] <- match(paste(newsplit[myperms[k,]],collapse=":"),colnames(ssamk$et))}
            pidx <- which(!is.na(newmatch))
            if(length(pidx)>0L){
              tnames[j] <- paste(newsplit[myperms[pidx,]],collapse=":")
            } else{stop("Cannot add new effects to formula when entering 'makessa' object.")}
          } else{stop("Cannot add new effects to formula when entering 'makessa' object.")}
        } # end for(j in idxna)
        colnames(et) <- tnames
      } # end if(any(is.na(mtidx)))
      oldnames <- rownames(ssamk$et)
      if(newnames[1]!=oldnames[1]){stop("Must include same response in formula when entering 'makessa' object.")}
      newnames <- newnames[2:mfdim[1]]
      oldnames <- oldnames[2:(ssamk$nxvar+1L)]
      
      # check order of predictors
      midx <- match(newnames,oldnames)
      if(any(is.na(midx))){stop("Cannot include new predictors in formula when entering 'makessa' object. \n Refit model using bigssa function.")}
      Etab <- matrix(0L,ssamk$nxvar,mfdim[2])
      Etab[midx,] <- et[-1,]
      rownames(Etab) <- oldnames
      colnames(Etab) <- tnames
      Etab <- rbind(0L,Etab)
      
    } # end if(formula==ssamk$formula)
    
    ### get info
    xnames <- oldnames
    tnames <- tnames
    xdim <- ssamk$xdim
    if(is.null(ssamk$gcvopts)){
      maxit <- 5
      gcvtol <- 10^-5
      alpha <- 1
    } else {
      maxit <- ssamk$gcvopts$maxit
      gcvtol <- ssamk$gcvopts$gcvtol
      alpha <- ssamk$gcvopts$alpha
    }
    if(is.null(ssamk$lambdas)){lambdas <- 10^-(9:0)} else {lambdas <- ssamk$lambdas}
    
    ### make design and penalty matrices
    dps <- ssadpm(ssamk$xvars,ssamk$type,ssamk$rks,ssamk$theknots,Etab)
    Knames <- colnames(dps$Kmat)
    jdim <- dim(dps$Jmats)
    Jnames <- colnames(dps$Jmat)[seq(1,jdim[2],by=ssamk$nknots)]
    
    ### get cross-product matrices
    wsqrt <- sqrt(ssamk$fweights*ssamk$weights)
    KtJ <- crossprod(dps$Kmat*(ssamk$fweights*ssamk$weights),dps$Jmats)
    KtK <- crossprod(dps$Kmat*wsqrt)
    JtJ <- crossprod(dps$Jmats*wsqrt)
    Kty <- crossprod(dps$Kmat*ssamk$weights,ssamk$xvars[[ssamk$nxvar+1]])
    Jty <- crossprod(dps$Jmats*ssamk$weights,ssamk$xvars[[ssamk$nxvar+1]])
    
    ### reml estimation of variance components
    if(!is.null(ssamk$random)){
      
      # some initializations
      XtZ <- t(makeZtX(ssamk$reinfo,dps$Kmat,dps$Jmats,ssamk$uidx,ssamk$ZtZy$Zmat))
      nbf <- length(Kty)
      cbf <- length(Jty)
      
      # iterative estimation (if remlalg!="none")
      if(ssamk$remlalg!="none"){
        if(jdim[2]>ssamk$nknots){
          gammat <- kronecker(rep(1,jdim[2]/ssamk$nknots),diag(ssamk$nknots))
          ktj <- KtJ%*%gammat
          if(length(ssamk$ZtZy$rdm)==1L){
            varhat <- remlri(ssamk$yty,rbind(Kty,crossprod(gammat,Jty)),matrix(ssamk$ZtZy[[1]]),rbind(cbind(KtK,ktj),cbind(t(ktj),crossprod(gammat,JtJ)%*%gammat)),
                             ssamk$ZtZy[[2]],rbind(XtZ[1:nbf,],crossprod(gammat,XtZ[(nbf+1):(nbf+cbf),])),ssamk$n[2]-ssamk$nknots,
                             tau=ssamk$remlopts$itau,imx=ssamk$remlopts$iter,tol=ssamk$remlopts$tol,alg=ssamk$remlalg)
          } else {
            varhat <- remlvc(ssamk$yty,rbind(Kty,crossprod(gammat,Jty)),matrix(ssamk$ZtZy[[1]]),rbind(cbind(KtK,ktj),cbind(t(ktj),crossprod(gammat,JtJ)%*%gammat)),
                             ssamk$ZtZy[[2]],rbind(XtZ[1:nbf,],crossprod(gammat,XtZ[(nbf+1):(nbf+cbf),])),ssamk$n[2]-ssamk$nknots,rdm=ssamk$ZtZy[[3]],
                             tau=ssamk$remlopts$itau,imx=ssamk$remlopts$iter,tol=ssamk$remlopts$tol,alg=ssamk$remlalg)
          }
          rm(ktj)
          
        } else {
          if(length(ssamk$ZtZy$rdm)==1L){
            varhat <- remlri(ssamk$yty,rbind(Kty,Jty),matrix(ssamk$ZtZy[[1]]),rbind(cbind(KtK,KtJ),cbind(t(KtJ),JtJ)),
                             ssamk$ZtZy[[2]],XtZ,ssamk$n[2]-ssamk$nknots,tau=ssamk$remlopts$itau,imx=ssamk$remlopts$iter,tol=ssamk$remlopts$tol,alg=ssamk$remlalg)
          } else{
            varhat <- remlvc(ssamk$yty,rbind(Kty,Jty),matrix(ssamk$ZtZy[[1]]),rbind(cbind(KtK,KtJ),cbind(t(KtJ),JtJ)),
                             ssamk$ZtZy[[2]],XtZ,ssamk$n[2]-ssamk$nknots,rdm=ssamk$ZtZy[[3]],tau=ssamk$remlopts$itau,imx=ssamk$remlopts$iter,tol=ssamk$remlopts$tol,alg=ssamk$remlalg)
          }
        }
        tau <- varhat$tau
        remlinfo <- varhat[3:5]
        if(!varhat$cnvg){warning("REML estimation algorithm failed to converge. \nTry a different algorithm (remlalg) or change REML estimation options (remlopts).")}
      } else {
        tau <- ssamk$remlopts$itau
        remlinfo <- NULL
      }
      nvc <- ssamk$ZtZy[[3]]
      names(tau) <- names(ssamk$ZtZy$rdm)
      
      # check for negative variance components
      if(any(tau<=0)){
        tau[tau<=0] <- 10^-4
        warning("REML estimation produced variance component estimates <= 0. \nOffending variance components are set to 10^-4 for estimation.")
      }
      
      # redefine crossproducts (part 1)
      if(is.matrix(ssamk$ZtZy[[2]])){
        # ZtZ is NOT diagonal (stored as matrix)
        csqrt <- pinvsm(ssamk$ZtZy[[2]]+diag(rep(1/tau,times=nvc)), iconst=-0.5)
        XtZc <- XtZ%*%csqrt
        if(nbf>1L) {
          Kty <- (Kty-(XtZc[1:nbf,])%*%(csqrt%*%ssamk$ZtZy[[1]]))
        } else {
          Kty <- (Kty-matrix(XtZc[1:nbf,],nbf,ncol(XtZc))%*%(csqrt%*%ssamk$ZtZy[[1]]))
        }
        Jty <- (Jty-(XtZc[(nbf+1):(nbf+cbf),])%*%(csqrt%*%ssamk$ZtZy[[1]]))
        ssamk$yty <- (ssamk$yty-crossprod(csqrt%*%ssamk$ZtZy[[1]]))
      } else {
        # ZtZ is diagonal (stored as vector)
        csqrt <- (ssamk$ZtZy[[2]]+rep(1/tau,times=nvc))^(-0.5)
        XtZc <- XtZ*rep(csqrt,each=(nbf+cbf))
        if(nbf>1L) {
          Kty <- (Kty-(XtZc[1:nbf,])%*%(ssamk$ZtZy[[1]]*csqrt))
        } else {
          Kty <- (Kty-matrix(XtZc[1:nbf,],nbf,ncol(XtZc))%*%(ssamk$ZtZy[[1]]*csqrt))
        }
        Jty <- (Jty-(XtZc[(nbf+1):(nbf+cbf),])%*%(ssamk$ZtZy[[1]]*csqrt))
        ssamk$yty <- (ssamk$yty-crossprod(ssamk$ZtZy[[1]]*csqrt))
      } # end if(is.matrix(ssamk$ZtZy[[2]]))
      
      # redefine crossproducts (part 2)
      if(nbf>1L){
        KtK <- (KtK-tcrossprod(XtZc[1:nbf,]))
        KtJ <- (KtJ-tcrossprod(XtZc[1:nbf,],XtZc[(nbf+1):(nbf+cbf),]))
      } else {
        nxtz <- ncol(XtZc)
        KtK <- (KtK-tcrossprod(matrix(XtZc[1:nbf,],nbf,nxtz)))
        KtJ <- (KtJ-tcrossprod(matrix(XtZc[1:nbf,],nbf,nxtz),XtZc[(nbf+1):(nbf+cbf),]))
      }
      JtJ <- (JtJ-tcrossprod(XtZc[(nbf+1):(nbf+cbf),]))
      
    } else {
      tau <- remlinfo <- NULL
    } # end if(!is.null(ssamk$random))
    
    
    ### initialize smoothing parameters
    nbf <- length(Kty)
    if(ssamk$nxvar>1L){
      if(is.null(ssamk$gammas[1])){
        gammas <- smartssa(dps$Qmats,Etab,Jnames,lambdas,Kty,Jty,KtK,KtJ,
                           JtJ,ssamk$nknots,ssamk$n[2],alpha,ssamk$yty,nbf)
      } else {
        gammas <- ssamk$gammas
      }
      if(length(gammas)<ssamk$nxvar){
        gidx <- which(rowSums(Etab)[2:(ssamk$nxvar+1L)]>0)
        mygammas <- rep(0,ssamk$nxvar)
        mygammas[gidx] <- gammas
        gammas <- mygammas
      }
    } else {
      gammas <- 1
    }
    
    ### estimate optimal smoothing parameters
    if(ssamk$skip.iter){
      cvg <- NA
      gamvec <- NULL
      for(j in 1:length(Jnames)){
        xi <- strsplit(Jnames[j],":")
        xidx <- match(xi[[1]],xnames)
        gamvec <- c(gamvec,prod(gammas[xidx]))
      }
      fxhat <- lamcoef(lambdas,gamvec,Kty,Jty,KtK,KtJ,JtJ,
                       dps$Qmats,ssamk$nknots,ssamk$n[2],alpha,ssamk$yty,nbf)
      fhat <- fxhat[[1]]
      dchat <- fhat[1:(nbf+ssamk$nknots)]
      yhat <- cbind(dps$Kmat,dps$Jmats%*%kronecker(gamvec,diag(ssamk$nknots)))%*%dchat
      sseval <- fhat[nbf+ssamk$nknots+1]
      effdf <- fhat[nbf+ssamk$nknots+2]
      if(sseval<=0){
        if(is.na(ssamk$rparm[1])){sseval <- sum((ssamk$xvars[[ssamk$nxvar+1]]-yhat)^2)} else{sseval <- .Machine$double.eps}
        warning("Approximated SSE is less than 0, so model 'info' may be invalid. \nTry reducing the number of knots or increasing lambda range (e.g., lambdas=10^-(5:0)).")
      }
      if(effdf>=(ssamk$n[2]-1L)){
        effdf <- ssamk$n[2]-.Machine$double.eps
        warning("Effective degrees of freedom exceeds n-1, so model 'info' may be invalid. \nTry reducing the number of knots or increasing lambda range (e.g., lambdas=10^-(5:0)).")
      }
      mevar <- sseval/(ssamk$n[2]-effdf)
      gcv <- ssamk$n[2]*sseval/((ssamk$n[2]-alpha*effdf)^2)
      newlam <- lambdas[fhat[nbf+ssamk$nknots+3]]
      aic <- ssamk$n[2]*(1+log(2*pi)) + ssamk$n[2]*log(sseval/ssamk$n[2]) + effdf*2
      bic <- ssamk$n[2]*(1+log(2*pi)) + ssamk$n[2]*log(sseval/ssamk$n[2]) + effdf*log(ssamk$n[2])
      csqrt <- sqrt(mevar)*fxhat[[2]]
      iter <- vtol <- NA
    } else {
      
      # iterate estimates of lambda and etas until convergence
      vtol <- 1
      gcv <- ssamk$yty
      iter <- 0L
      cvg <- FALSE
      while(vtol>gcvtol && iter<maxit && min(c(gammas,gcv))>0) {

        if(ssamk$nxvar==1L){
          nqmat <- matrix(0,nbf+ssamk$nknots,nbf+ssamk$nknots)
          nqmat[(nbf+1):(ssamk$nknots+nbf),(nbf+1):(ssamk$nknots+nbf)] <- ssamk$n[2]*dps$Qmats
          if(length(lambdas)>1){
            newlam <- lamloop(lambdas,gammas,Kty,Jty,KtK,KtJ,JtJ,dps$Qmats,
                              ssamk$nknots,ssamk$n[2],alpha,ssamk$yty,nbf)
          } else{
            newlam <- lambdas
          }
          gcvopt <- nlm(f=gcvcss,p=log(newlam),yty=ssamk$yty,xtx=rbind(cbind(KtK,KtJ),
                        cbind(t(KtJ),JtJ)),xty=rbind(Kty,Jty),nqmat=nqmat,
                        ndpts=ssamk$n[2],alpha=alpha)
          gcv <- gcvopt$min
          newlam <- exp(gcvopt$est)
          vtol <- 0
        } else{
          
          # step 1: find optimal lambda given gammas
          gamvec <- NULL
          for(j in 1:length(Jnames)){
            xi <- strsplit(Jnames[j],":")
            xidx <- match(xi[[1]],xnames)
            gamvec <- c(gamvec,prod(gammas[xidx]))
          }
          newlam <- lamloop(lambdas,gamvec,Kty,Jty,KtK,KtJ,JtJ,dps$Qmats,
                            ssamk$nknots,ssamk$n[2],alpha,ssamk$yty,nbf)
          
          # step 2: find optimal etas given lambda
          gcvopt <- nlm(f=gcvssa,p=log(gammas),yty=ssamk$yty,KtK=KtK,KtJ=KtJ,Kty=Kty,
                        Jty=Jty,JtJ=JtJ,Qmats=dps$Qmats,ndpts=ssamk$n[2],alpha=alpha,
                        nknots=ssamk$nknots,newlam=newlam,nbf=nbf,xnames=xnames,Jnames=Jnames)
          newgcv <- gcvopt$min
          
          # step 3: check for convergence
          gammas <- exp(gcvopt$est)
          vtol <- (gcv-newgcv)/gcv
          gcv <- newgcv
          iter <- iter+1
          
        } # end if(ncol(dps$Qmats)==ssamk$nknots)
      } # end while(vtol>gcvtol && iter<maxit && min(c(gammas,gcv))>0)
      
      # get final estimates
      gamvec <- NULL
      for(j in 1:length(Jnames)){
        xi <- strsplit(Jnames[j],":")
        xidx <- match(xi[[1]],xnames)
        gamvec <- c(gamvec,prod(gammas[xidx]))
      }
      fxhat <- lamcoef(newlam,gamvec,Kty,Jty,KtK,KtJ,JtJ,dps$Qmats,
                       ssamk$nknots,ssamk$n[2],alpha,ssamk$yty,nbf)
      fhat <- fxhat[[1]]
      dchat <- fhat[1:(nbf+ssamk$nknots)]
      yhat <- cbind(dps$Kmat,dps$Jmats%*%kronecker(gamvec,diag(ssamk$nknots)))%*%dchat
      sseval <- fhat[nbf+ssamk$nknots+1]
      effdf <- fhat[nbf+ssamk$nknots+2]
      if(sseval<=0){
        if(is.na(ssamk$rparm[1])){sseval <- sum((ssamk$xvars[[ssamk$nxvar+1]]-yhat)^2)} else{sseval <- .Machine$double.eps}
        warning("Approximated SSE is less than 0, so model 'info' may be invalid. \nTry reducing the number of knots or increasing lambda range (e.g., lambdas=10^-(5:0)).")
      }
      if(effdf>=(ssamk$n[2]-1L)){
        effdf <- ssamk$n[2]-.Machine$double.eps
        warning("Effective degrees of freedom exceeds n-1, so model 'info' may be invalid. \nTry reducing the number of knots or increasing lambda range (e.g., lambdas=10^-(5:0)).")
      }
      mevar <- sseval/(ssamk$n[2]-effdf)
      gcv <- ssamk$n[2]*sseval/((ssamk$n[2]-alpha*effdf)^2)
      aic <- ssamk$n[2]*(1+log(2*pi)) + ssamk$n[2]*log(sseval/ssamk$n[2]) + effdf*2
      bic <- ssamk$n[2]*(1+log(2*pi)) + ssamk$n[2]*log(sseval/ssamk$n[2]) + effdf*log(ssamk$n[2])
      csqrt <- sqrt(mevar)*fxhat[[2]]
      
    } # end if(ssamk$skip.iter)
    
    ### posterior variance
    pse=NA
    if(ssamk$se.fit){pse <- sqrt(postvar(dps$Kmat,dps$Jmats%*%kronecker(gamvec,diag(ssamk$nknots)),csqrt))}
    
    ### calculate vaf
    mval <- ssamk$ysm/ssamk$n[2]
    if(!is.null(ssamk$random)){
      if(is.matrix(ssamk$ZtZy[[2]])){
        mval <- NA
      } else{
        mval <- sum(sqrt(1-(ssamk$ZtZy[[2]]/(ssamk$ZtZy[[2]]+rep(1/tau,times=nvc))))*ssamk$ZtZy[[1]])/ssamk$n[2]
      }
    }
    vaf <- 1-sseval/(ssamk$yty-ssamk$n[2]*(mval^2))
    
    ### retransform predictors
    for(k in 1:ssamk$nxvar){
      if(any(ssamk$type[[k]]==c("cub","cub0","per"))){
        ssamk$xvars[[k]] <- as.matrix(ssamk$xvars[[k]]*(ssamk$xrng[[k]][2]-ssamk$xrng[[k]][1])+ssamk$xrng[[k]][1])
      } else if(ssamk$type[[k]]=="ord"){
        ssamk$xvars[[k]] <- as.matrix(ssamk$flvls[[k]][ssamk$xvars[[k]]])
        if(!is.na(ssamk$rparm[1])){ssamk$xorig[[k]] <- as.matrix(ssamk$flvls[[k]][ssamk$xorig[[k]]])}
      } else if(ssamk$type[[k]]=="nom"){
        ssamk$xvars[[k]] <- as.matrix(ssamk$flvls[[k]][ssamk$xvars[[k]]])
        if(!is.na(ssamk$rparm[1])){ssamk$xorig[[k]] <- as.matrix(ssamk$flvls[[k]][ssamk$xorig[[k]]])}
      }
    } # end for(k in 1:ssamk$nxvar)
    if(!is.na(ssamk$rparm[1])){
      yunique <- ssamk$xvars[[ssamk$nxvar+1]]/as.numeric(ssamk$fweight); yvar=ssamk$yorig
      xunique <- ssamk$xvars[1:ssamk$nxvar]
      ssamk$xvars <- ssamk$xorig
    } else{
      xunique <- yunique <- NA
      yvar <- ssamk$xvars[[ssamk$nxvar+1]]
      ssamk$xvars <- ssamk$xvars[1:ssamk$nxvar]
    }
    
    ### check for random (get blup)
    bhat <- NULL
    if(!is.null(ssamk$random)){
      if(jdim[2]>ssamk$nknots){
        gamvec <- NULL
        for(j in 1:length(Jnames)){
          xi <- strsplit(Jnames[j],":")
          xidx <- match(xi[[1]],xnames)
          gamvec <- c(gamvec,prod(gammas[xidx]))
        }
        gammat <- kronecker(gamvec,diag(ssamk$nknots))
        XtZ <- rbind(XtZ[1:nbf,],crossprod(gammat,XtZ[(nbf+1):(nbf+cbf),]))
      }
      bhat <- ssblup(ssamk$ZtZy[[1]],t(XtZ),ssamk$ZtZy[[2]],dchat,ssamk$ZtZy[[3]],tau)
    }
    
    ### collect results
    names(gammas) <- xnames
    if(ssamk$skip.iter==FALSE && vtol<=gcvtol){cvg <- TRUE}
    ndf <- data.frame(n=ssamk$n[2],df=effdf,row.names="")
    modelspec <- list(myknots=ssamk$theknots,rparm=ssamk$rparm,lambda=newlam,
                      gammas=gammas,gcvopts=ssamk$gcvopts,nxvar=xdim,xrng=ssamk$xrng,
                      flvls=ssamk$flvls,tpsinfo=ssamk$tpsinfo,iter=iter,vtol=vtol,
                      coef=dchat,coef.csqrt=csqrt,Etab=Etab,Knames=Knames,
                      fweights=ssamk$fweights,weights=ssamk$weights,remlalg=ssamk$remlalg,
                      remlopts=ssamk$remlopts,remlinfo=remlinfo)
    ssafit <- list(fitted.values=yhat,se.fit=pse,yvar=yvar,xvars=ssamk$xvars,
                   type=ssamk$type,yunique=yunique,xunique=xunique,sigma=sqrt(mevar),
                   ndf=ndf,info=c(gcv=gcv,rsq=vaf,aic=aic,bic=bic),modelspec=modelspec,
                   converged=cvg,tnames=tnames,random=ssamk$random,tau=tau,blup=bhat)
    return(ssafit)
    
  }