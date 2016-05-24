ssgwork <-
  function(formula,ssgmk){
    
    ### check for negbin with missing dispersion
    if(ssgmk$family=="negbin"){
      if(is.null(ssgmk$dispersion)){
        ssgmk$family <- "poisson"
        se.lpold <- ssgmk$se.lp
        ssgmk$se.lp <- TRUE
        pmod <- ssgwork(formula,ssgmk)
        top <- sum(pmod$fitted*(1-pmod$fitted*(pmod$se.lp^2))*pmod$modelspec$fweights)
        bot <- sum(ssgmk$yty/pmod$fitted)-2*ssgmk$ysm+sum(pmod$fitted*pmod$modelspec$fweights)-(pmod$ndf[1]-pmod$ndf[2])
        ssgmk$dispersion <- as.numeric(bot/top)
        rm(pmod)
        ssgmk$family <- "negbin"
        ssgmk$se.lp <- se.lpold
        estdisp <- TRUE
      } else {
        estdisp <- FALSE
      }
    } else {
      estdisp <- FALSE
    }
    
    ### check inputs
    if(formula==ssgmk$formula){
      Etab <- ssgmk$et
      mfdim <- dim(Etab)
      tnames <- colnames(Etab)
      oldnames <- rownames(ssgmk$et)[2:(ssgmk$nxvar+1L)]
    } else{
      
      # get formula
      trmfrm <- terms.formula(formula)
      ychk <- attr(trmfrm,"response")
      if(ychk!=1L){stop("Must include same response in formula when entering 'makessa' object.")}
      et <- attr(trmfrm,"factors")
      mfdim <- dim(et)
      newnames <- rownames(et)
      tnames <- colnames(et)
      mtidx <- match(tnames,colnames(ssgmk$et))
      if(any(is.na(mtidx))){
        idxna <- which(is.na(mtidx))
        for(j in idxna){
          newsplit <- unlist(strsplit(tnames[j],":"))
          lnwsp <- length(newsplit)
          if(lnwsp==1L){
            stop("Cannot add new effects to formula when entering 'makessa' object.")
          } else if(lnwsp==2L){
            newmatch <- match(paste(newsplit[2],newsplit[1],sep=":"),colnames(ssgmk$et))
            if(is.na(newmatch)){
              stop("Cannot add new effects to formula when entering 'makessa' object.")
            } else{
              tnames[j] <- paste(newsplit[2],newsplit[1],sep=":")
            }
          } else if(lnwsp==3L){
            myperms <- matrix(c(3,2,3,1,2,2,3,1,3,1,1,1,2,2,3),ncol=3)
            newmatch <- rep(NA,5)
            for(k in 1:5){newmatch[k] <- match(paste(newsplit[myperms[k,]],collapse=":"),colnames(ssgmk$et))}
            pidx <- which(is.na(newmatch)==FALSE)
            if(length(pidx)>0L){
              tnames[j] <- paste(newsplit[myperms[pidx,]],collapse=":")
            } else{stop("Cannot add new effects to formula when entering 'makessa' object.")}
          } else{stop("Cannot add new effects to formula when entering 'makessa' object.")}
        } # end for(j in idxna)
        colnames(et) <- tnames
      } # end if(any(is.na(mtidx)))
      oldnames <- rownames(ssgmk$et)
      if(newnames[1]!=oldnames[1]){stop("Must include same response in formula when entering 'makessa' object.")}
      newnames <- newnames[2:mfdim[1]]
      oldnames <- oldnames[2:(ssgmk$nxvar+1L)]
      
      # check order of predictors
      midx <- match(newnames,oldnames)
      if(any(is.na(midx))){stop("Cannot include new predictors in formula when entering 'makessa' object. \n Refit model using bigssa function.")}
      Etab <- matrix(0L,ssgmk$nxvar,mfdim[2])
      Etab[midx,] <- et[-1,]
      rownames(Etab) <- oldnames
      colnames(Etab) <- tnames
      Etab <- rbind(0L,Etab)
      
    } # end if(formula==ssgmk$formula)
    
    ### get info
    xnames <- oldnames
    tnames <- tnames
    xdim <- ssgmk$xdim
    if(is.null(ssgmk$gcvopts)){
      maxit <- 5
      gcvtol <- 10^-5
      alpha <- 1
      inmaxit <- 100
      intol <- 10^-5
      insub <- 10000
    } else {
      maxit <- ssgmk$gcvopts$maxit
      gcvtol <- ssgmk$gcvopts$gcvtol
      alpha <- ssgmk$gcvopts$alpha
      inmaxit <- ssgmk$gcvopts$inmaxit
      intol <- ssgmk$gcvopts$intol
      insub <- ssgmk$gcvopts$insub
    }
    if(is.null(ssgmk$lambdas)){lambdas <- 10^-(9:0)} else {lambdas <- ssgmk$lambdas}
    
    ### make design and penalty matrices
    dps <- ssadpm(ssgmk$xvars,ssgmk$type,ssgmk$rks,ssgmk$theknots,Etab)
    Knames <- colnames(dps$Kmat)
    jdim <- dim(dps$Jmats)
    Jnames <- colnames(dps$Jmat)[seq(1,jdim[2],by=ssgmk$nknots)]
    
    ### initialize smoothing parameters
    nbf <- ncol(dps$Kmat)
    if(ssgmk$nxvar>1L){
      if(is.null(ssgmk$gammas[1])){
        gammas <- smartssg(Etab,lambdas,ssgmk$family,dps$Kmat,dps$Jmats,Jnames,
                           ssgmk$xvars[[ssgmk$nxvar+1]],dps$Qmats,ssgmk$nknots,
                           ssgmk$n[2],alpha,ssgmk$yty,nbf,ssgmk$fweights,ssgmk$weights,
                           inmaxit,intol,insub,ssgmk$dispersion,ssgmk$gcvtype)
      } else {
        gammas <- ssgmk$gammas
      }
      if(length(gammas)<ssgmk$nxvar){
        gidx <- which(rowSums(Etab)[2:(ssgmk$nxvar+1L)]>0)
        mygammas <- rep(0,ssgmk$nxvar)
        mygammas[gidx] <- gammas
        gammas <- mygammas
      }
    } else {
      gammas <- 1
    }
    
    ### estimate optimal smoothing parameters
    if(ssgmk$skip.iter){
      cvg <- NA
      gamvec <- NULL
      for(j in 1:length(Jnames)){
        xi <- strsplit(Jnames[j],":")
        xidx <- match(xi[[1]],xnames)
        gamvec <- c(gamvec,prod(gammas[xidx]))
      }
      fxhat <- lamcoefg(lambdas,gamvec,ssgmk$family,dps$Kmat,dps$Jmats,
                        ssgmk$xvars[[ssgmk$nxvar+1]],dps$Qmats,ssgmk$nknots,
                        ssgmk$n[2],alpha,ssgmk$yty,nbf,ssgmk$fweights,ssgmk$weights,
                        inmaxit,intol,insub,ssgmk$dispersion,ssgmk$gcvtype)
      fhat <- fxhat[[1]]
      dchat <- fhat[1:(nbf+ssgmk$nknots)]
      gcv <- fhat[nbf+ssgmk$nknots+1]
      effdf <- fhat[nbf+ssgmk$nknots+2]
      newlam <- lambdas[fhat[nbf+ssgmk$nknots+3]]
      yhat <- cbind(dps$Kmat,dps$Jmats%*%kronecker(gamvec,diag(ssgmk$nknots)))%*%dchat
      if(ssgmk$family=="binomial"){
        mevar <- 1
        fitvals <- 1/(1+exp(-yhat))
        n2LL <- (-2)*(crossprod(ssgmk$xvars[[ssgmk$nxvar+1]],yhat) + sum(ssgmk$weights*ssgmk$fweights*log(1/(1+exp(yhat)))) + ssgmk$slogy)
      } else if(ssgmk$family=="poisson"){
        mevar <- 1
        fitvals <- exp(yhat)
        n2LL <- (-2)*(crossprod(ssgmk$xvars[[ssgmk$nxvar+1]],yhat) - sum(fitvals*ssgmk$fweights) - ssgmk$slogy)
      } else if(ssgmk$family=="Gamma"){
        mevar <- fhat[nbf+ssgmk$nknots+4]/(ssgmk$n[2]-effdf)
        fitvals <- 1/yhat
        n2LL <- (-2)*( (1/mevar)*(crossprod(ssgmk$xvars[[ssgmk$nxvar+1]],-yhat) - sum(log(fitvals)*ssgmk$fweights)) + ((1/mevar)-1)*ssgmk$slogy - ssgmk$n[2]*(1/mevar)*log(mevar) - ssgmk$n[2]*log(gamma(1/mevar))  )
      } else if(ssgmk$family=="inverse.gaussian"){
        mevar <- fhat[nbf+ssgmk$nknots+4]/(ssgmk$n[2]-effdf)
        fitvals <- sqrt(1/yhat)
        n2LL <- (-2)*(-(1/2)*( (1/mevar)*((-2)*crossprod(ssgmk$xvars[[ssgmk$nxvar+1]],-yhat/2) - 2*sum(ssgmk$fweights/fitvals) + ssgmk$slogy[1]) + ssgmk$n[2]*log(mevar) + ssgmk$slogy[2] + ssgmk$n[2]*log(2*pi) ))
      } else if(ssgmk$family=="negbin"){
        mevar <- 1
        fitvals <- exp(yhat)
        if(is.na(ssgmk$rparm[1])){
          if(estdisp){size <- nbmle(ssgmk$xvars[[ssgmk$nxvar+1]],fitvals,fhat[nbf+ssgmk$nknots+4],ssgmk$fweights,ssgmk$n[2],ssgmk$xvars[[ssgmk$nxvar+1]])} else {size <- fhat[nbf+ssgmk$nknots+4]}
          slogy <- sum(lgamma(ssgmk$xvars[[ssgmk$nxvar+1]]+size))
        } else{
          if(estdisp){size <- nbmle(ssgmk$yorig,fitvals,fhat[nbf+ssgmk$nknots+4],ssgmk$fweights,ssgmk$n[2],ssgmk$xvars[[ssgmk$nxvar+1]])} else {size <- fhat[nbf+ssgmk$nknots+4]}
          slogy <- sum(lgamma(ssgmk$yorig+size))
        }
        n2LL <- (-2)*( crossprod(ssgmk$xvars[[ssgmk$nxvar+1]],yhat) - ssgmk$ysm*log(size) - sum((ssgmk$xvars[[ssgmk$nxvar+1]] + size*ssgmk$fweights)*log(1+fitvals/size)) + slogy - ssgmk$slogy - ssgmk$n[2]*lgamma(size) )
      }
      if(effdf>=(ssgmk$n[2]-1L)){
        effdf <- ssgmk$n[2]-.Machine$double.eps
        warning("Effective degrees of freedom exceeds n-1, so model 'info' may be invalid. \nTry reducing the number of knots or increasing lambda range (e.g., lambdas=10^-(5:0)).")
      }
      aic <- n2LL+effdf*2
      bic <- n2LL+effdf*log(ssgmk$n[2])
      csqrt <- sqrt(mevar)*fxhat[[2]]
      if(ssgmk$family=="negbin"){mevar <- 1/size}
      iter <- vtol <- NA
    } else {
      
      # iterate estimates of lambda and etas until convergence
      vtol <- 1
      gcv <- sum(ssgmk$yty)
      iter <- 0L
      cvg <- FALSE
      while(vtol>gcvtol && iter<maxit && min(gammas)>0) {
        
        if(ssgmk$nxvar==1L){
          if(length(lambdas)>1){
            newlam <- lamloopg(lambdas,gammas,ssgmk$family,dps$Kmat,dps$Jmats,
                               ssgmk$xvars[[ssgmk$nxvar+1]],dps$Qmats,ssgmk$nknots,
                               ssgmk$n[2],alpha,ssgmk$yty,nbf,ssgmk$fweights,ssgmk$weights,
                               inmaxit,intol,insub,ssgmk$dispersion,ssgmk$gcvtype)
          } else{
            newlam <- lambdas
          }
          gcvopt <- nlm(f=gcvgss,p=log(newlam),family=ssgmk$family,Kmat=dps$Kmat,Jmat=dps$Jmats,
                        yvar=ssgmk$xvars[[ssgmk$nxvar+1]],Qmat=dps$Qmats,nknots=ssgmk$nknots,
                        ndpts=ssgmk$n[2],alpha=alpha,yty=ssgmk$yty,nbf=nbf,
                        fweights=ssgmk$fweights,weights=ssgmk$weights,maxit=inmaxit,intol=intol,
                        subsamp=insub,dispersion=ssgmk$dispersion,gcvtype=ssgmk$gcvtype)
          gcv <- gcvopt$min
          oldlam <- newlam
          newlam <- exp(gcvopt$est)
          if(ssgmk$family=="negbin" & estdisp){
            fhat <- lamcoefg(newlam,1,ssgmk$family,dps$Kmat,dps$Jmats,
                             ssgmk$xvars[[ssgmk$nxvar+1]],dps$Qmats,ssgmk$nknots,
                             ssgmk$n[2],alpha,ssgmk$yty,nbf,ssgmk$fweights,ssgmk$weights,
                             inmaxit,intol,insub,ssgmk$dispersion,ssgmk$gcvtype)[[1]]
            mu0 <- exp(cbind(dps$Kmat,dps$Jmats%*%kronecker(1,diag(ssgmk$nknots)))%*%fhat[1:(nbf+ssgmk$nknots)])
            size <- nbmle(ssgmk$yorig,mu0,1/ssgmk$dispersion,ssgmk$fweights,ssgmk$n[2],ssgmk$xvars[[ssgmk$nxvar+1]])
            vtol <- abs((oldlam-newlam)/oldlam)
            ssgmk$dispersion <- 1/size
            lambdas <- newlam
            oldlam <- newlam
            iter <- iter+1L
          } else {
            vtol <- 0
          }
        } else{
          
          # step 1: find optimal lambda given gammas
          gamvec <- NULL
          for(j in 1:length(Jnames)){
            xi <- strsplit(Jnames[j],":")
            xidx <- match(xi[[1]],xnames)
            gamvec <- c(gamvec,prod(gammas[xidx]))
          }
          newlam <- lamloopg(lambdas,gamvec,ssgmk$family,dps$Kmat,dps$Jmats,
                             ssgmk$xvars[[ssgmk$nxvar+1]],dps$Qmats,ssgmk$nknots,
                             ssgmk$n[2],alpha,ssgmk$yty,nbf,ssgmk$fweights,ssgmk$weights,
                             inmaxit,intol,insub,ssgmk$dispersion,ssgmk$gcvtype)
          
          # step 2: find optimal etas given lambda
          gcvopt <- nlm(f=gcvssg,p=log(gammas),newlam=newlam,family=ssgmk$family,
                        Kmat=dps$Kmat,Jmats=dps$Jmats,yvar=ssgmk$xvars[[ssgmk$nxvar+1]],
                        Qmats=dps$Qmats,nknots=ssgmk$nknots,ndpts=ssgmk$n[2],alpha=alpha,
                        yty=ssgmk$yty,nbf=nbf,fweights=ssgmk$fweights,weights=ssgmk$weights,
                        xnames=xnames,Jnames=Jnames,maxit=inmaxit,intol=intol,subsamp=insub,
                        dispersion=ssgmk$dispersion,gcvtype=ssgmk$gcvtype)
          newgcv <- gcvopt$min
          
          # step 3: check for convergence
          gammas <- exp(gcvopt$est)
          vtol <- (gcv-newgcv)/gcv
          gcv <- newgcv
          iter <- iter+1
          
          # step 4: update dispersion
          if(ssgmk$family=="negbin" & estdisp){
            gamvec <- NULL
            for(j in 1:length(Jnames)){
              xi <- strsplit(Jnames[j],":")
              xidx <- match(xi[[1]],xnames)
              gamvec <- c(gamvec,prod(gammas[xidx]))
            }
            fhat <- lamcoefg(newlam,gamvec,ssgmk$family,dps$Kmat,dps$Jmats,
                             ssgmk$xvars[[ssgmk$nxvar+1]],dps$Qmats,ssgmk$nknots,
                             ssgmk$n[2],alpha,ssgmk$yty,nbf,ssgmk$fweights,ssgmk$weights,
                             inmaxit,intol,insub,ssgmk$dispersion,ssgmk$gcvtype)[[1]]
            mu0 <- exp(cbind(dps$Kmat,dps$Jmats%*%kronecker(gamvec,diag(ssgmk$nknots)))%*%fhat[1:(nbf+ssgmk$nknots)])
            size <- nbmle(ssgmk$yorig,mu0,1/ssgmk$dispersion,ssgmk$fweights,ssgmk$n[2],ssgmk$xvars[[ssgmk$nxvar+1]])
            ssgmk$dispersion <- 1/size
          }
          
        } # end if(ncol(dps$Qmats)==ssgmk$nknots)
      } # end while(vtol>gcvtol && iter<maxit && min(c(gammas,gcv))>0)
      
      # get final estimates
      gamvec <- NULL
      for(j in 1:length(Jnames)){
        xi <- strsplit(Jnames[j],":")
        xidx <- match(xi[[1]],xnames)
        gamvec <- c(gamvec,prod(gammas[xidx]))
      }
      fxhat <- lamcoefg(newlam,gamvec,ssgmk$family,dps$Kmat,dps$Jmats,
                        ssgmk$xvars[[ssgmk$nxvar+1]],dps$Qmats,ssgmk$nknots,
                        ssgmk$n[2],alpha,ssgmk$yty,nbf,ssgmk$fweights,ssgmk$weights,
                        inmaxit,intol,insub,ssgmk$dispersion,ssgmk$gcvtype)
      fhat <- fxhat[[1]]
      dchat <- fhat[1:(nbf+ssgmk$nknots)]
      effdf <- fhat[nbf+ssgmk$nknots+2]
      yhat <- cbind(dps$Kmat,dps$Jmats%*%kronecker(gamvec,diag(ssgmk$nknots)))%*%dchat
      if(ssgmk$family=="binomial"){
        mevar <- 1
        fitvals <- 1/(1+exp(-yhat))
        n2LL <- (-2)*(crossprod(ssgmk$xvars[[ssgmk$nxvar+1]],yhat) + sum(ssgmk$weights*ssgmk$fweights*log(1/(1+exp(yhat)))) + ssgmk$slogy)
      } else if(ssgmk$family=="poisson"){
        mevar <- 1
        fitvals <- exp(yhat)
        n2LL <- (-2)*(crossprod(ssgmk$xvars[[ssgmk$nxvar+1]],yhat) - sum(fitvals*ssgmk$fweights) - ssgmk$slogy)
      } else if(ssgmk$family=="Gamma"){
        mevar <- fhat[nbf+ssgmk$nknots+4]/(ssgmk$n[2]-effdf)
        fitvals <- 1/yhat
        n2LL <- (-2)*( (1/mevar)*(crossprod(ssgmk$xvars[[ssgmk$nxvar+1]],-yhat) - sum(log(fitvals)*ssgmk$fweights)) + ((1/mevar)-1)*ssgmk$slogy - ssgmk$n[2]*(1/mevar)*log(mevar) - ssgmk$n[2]*log(gamma(1/mevar))  )
      } else if(ssgmk$family=="inverse.gaussian"){
        mevar <- fhat[nbf+ssgmk$nknots+4]/(ssgmk$n[2]-effdf)
        fitvals <- sqrt(1/yhat)
        n2LL <- (-2)*(-(1/2)*( (1/mevar)*((-2)*crossprod(ssgmk$xvars[[ssgmk$nxvar+1]],-yhat/2) - 2*sum(ssgmk$fweights/fitvals) + ssgmk$slogy[1]) + ssgmk$n[2]*log(mevar) + ssgmk$slogy[2] + ssgmk$n[2]*log(2*pi) ))
      } else if(ssgmk$family=="negbin"){
        mevar <- 1
        fitvals <- exp(yhat)
        if(is.na(ssgmk$rparm[1])){
          if(estdisp){size <- nbmle(ssgmk$xvars[[ssgmk$nxvar+1]],fitvals,fhat[nbf+ssgmk$nknots+4],ssgmk$fweights,ssgmk$n[2],ssgmk$xvars[[ssgmk$nxvar+1]])} else {size <- fhat[nbf+ssgmk$nknots+4]}
          slogy <- sum(lgamma(ssgmk$xvars[[ssgmk$nxvar+1]]+size))
        } else{
          if(estdisp){size <- nbmle(ssgmk$yorig,fitvals,fhat[nbf+ssgmk$nknots+4],ssgmk$fweights,ssgmk$n[2],ssgmk$xvars[[ssgmk$nxvar+1]])} else {size <- fhat[nbf+ssgmk$nknots+4]}
          slogy <- sum(lgamma(ssgmk$yorig+size))
        }
        n2LL <- (-2)*( crossprod(ssgmk$xvars[[ssgmk$nxvar+1]],yhat) - ssgmk$ysm*log(size) - sum((ssgmk$xvars[[ssgmk$nxvar+1]]+size*ssgmk$fweights)*log(1+fitvals/size)) + slogy - ssgmk$slogy - ssgmk$n[2]*lgamma(size) )
      }
      if(effdf>=(ssgmk$n[2]-1L)){
        effdf <- ssgmk$n[2]-.Machine$double.eps
        warning("Effective degrees of freedom exceeds n-1, so model 'info' may be invalid. \nTry reducing the number of knots or increasing lambda range (e.g., lambdas=10^-(5:0)).")
      }
      aic <- n2LL+effdf*2
      bic <- n2LL+effdf*log(ssgmk$n[2])
      csqrt <- sqrt(mevar)*fxhat[[2]]
      if(ssgmk$family=="negbin"){mevar <- 1/size}
      
    } # end if(ssgmk$skip.iter)
    
    ### posterior variance
    pse <- NA
    if(ssgmk$se.lp){pse <- sqrt(postvar(dps$Kmat,dps$Jmats%*%kronecker(gamvec,diag(ssgmk$nknots)),csqrt))}
    
    ### calculate vaf
    mval <- ssgmk$ysm/ssgmk$n[2]
    if(ssgmk$family=="binomial"){
      vaf <- 1-(sum(ssgmk$yty)-2*sum(ssgmk$xvars[[ssgmk$nxvar+1]]*fitvals/ssgmk$weights)+crossprod(fitvals^2,ssgmk$fweights))/(sum(ssgmk$yty)-ssgmk$n[2]*(mval^2))
    } else{
      vaf <- 1-(sum(ssgmk$yty)-2*sum(ssgmk$xvars[[ssgmk$nxvar+1]]*fitvals)+crossprod(fitvals^2,ssgmk$fweights))/(sum(ssgmk$yty)-ssgmk$n[2]*(mval^2))
    }
    
    ### retransform predictors
    for(k in 1:ssgmk$nxvar){
      if(any(ssgmk$type[[k]]==c("cub","cub0","per"))){
        ssgmk$xvars[[k]] <- as.matrix(ssgmk$xvars[[k]]*(ssgmk$xrng[[k]][2]-ssgmk$xrng[[k]][1])+ssgmk$xrng[[k]][1])
      } else if(ssgmk$type[[k]]=="ord"){
        ssgmk$xvars[[k]] <- as.matrix(ssgmk$flvls[[k]][ssgmk$xvars[[k]]])
        if(!is.na(ssgmk$rparm[1])){ssgmk$xorig[[k]] <- as.matrix(ssgmk$flvls[[k]][ssgmk$xorig[[k]]])}
      } else if(ssgmk$type[[k]]=="nom"){
        ssgmk$xvars[[k]] <- as.matrix(ssgmk$flvls[[k]][ssgmk$xvars[[k]]])
        if(!is.na(ssgmk$rparm[1])){ssgmk$xorig[[k]] <- as.matrix(ssgmk$flvls[[k]][ssgmk$xorig[[k]]])}
      }
    } # end for(k in 1:ssgmk$nxvar)
    if(!is.na(ssgmk$rparm[1])){
      yunique <- ssgmk$xvars[[ssgmk$nxvar+1]]/as.numeric(ssgmk$fweight)
      yvar <- ssgmk$yorig
      xunique <- ssgmk$xvars[1:ssgmk$nxvar]
      ssgmk$xvars <- ssgmk$xorig
      ssgmk$weights <- ssgmk$worig
    } else{
      xunique <- yunique <- NA
      yvar <- ssgmk$xvars[[ssgmk$nxvar+1]]
      ssgmk$xvars <- ssgmk$xvars[1:ssgmk$nxvar]
    }
    
    ### collect results
    names(gammas) <- xnames
    if(ssgmk$skip.iter==FALSE && vtol<=gcvtol){cvg <- TRUE}
    ndf <- data.frame(n=ssgmk$n[2],df=effdf,row.names="")
    modelspec <- list(myknots=ssgmk$theknots,rparm=ssgmk$rparm,lambda=newlam,
                      gammas=gammas,gcvopts=ssgmk$gcvopts,nxvar=xdim,xrng=ssgmk$xrng,
                      flvls=ssgmk$flvls,tpsinfo=ssgmk$tpsinfo,iter=iter,vtol=vtol,
                      coef=dchat,coef.csqrt=csqrt,Etab=Etab,Knames=Knames,
                      fweights=ssgmk$fweights,weights=ssgmk$weights,gcvtype=ssgmk$gcvtype)
    ssgfit <- list(fitted.values=fitvals,linear.predictors=yhat,se.lp=pse,yvar=yvar,
                   xvars=ssgmk$xvars,type=ssgmk$type,yunique=yunique,xunique=xunique,dispersion=mevar,
                   ndf=ndf,info=c(gcv=gcv,rsq=vaf,aic=aic,bic=bic),modelspec=modelspec,
                   converged=cvg,tnames=tnames,family=ssgmk$family)
    return(ssgfit)
    
  }