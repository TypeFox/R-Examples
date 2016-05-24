predict.bigssg <-
  function(object,newdata=NULL,se.lp=FALSE,include=object$tnames,
           effect=c("all","0","lin","non"),includeint=FALSE,
           design=FALSE,smoothMatrix=FALSE,...) {
    ###### Predicts for class "bigssg" objects
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: March 6, 2016
    
    ### check newdata
    effect <- effect[1]
    if(!any(effect==c("all","0","lin","non"))){stop("Must set 'effect' to one of four specified options.")}
    if(effect=="0"){
      yhat <- object$modelspec$coef[1]
      if(object$family=="binomial"){
        fitvals <- 1/(1+exp(-yhat))
      } else if(object$family=="poisson"){
        fitvals <- exp(yhat)
      } else if(object$family=="Gamma"){
        fitvals <- 1/yhat
      } else if(object$family=="inverse.gaussian"){
        fitvals <- sqrt(1/yhat)
      } else if(object$family=="negbin"){
        fitvals <- exp(yhat)
      }
      pred <- list(fitted.values=as.numeric(fitvals),linear.predictors=as.numeric(yhat))
      if(se.lp){
        pse <- sqrt(sum(object$modelspec$coef.csqrt[1,]^2))
        pred <- c(pred,list(se.lp=pse))
      } 
      return(pred)
    }
    xnames <- names(object$xvars)
    nxvar <- length(object$xvars)
    if(is.null(newdata)){
      newdata <- object$xvars
    } else {
      newdata <- as.list(newdata)
      newDnames <- names(newdata)
      if(is.null(newDnames)){stop("Variable names in 'newdata' must match those of 'xvars' element of input 'object'.")}
      midx <- match(newDnames,xnames)
      widx <- which(is.na(midx))
      if(length(widx)>0){
        if(length(widx)==length(newdata)){stop("Variable names in 'newdata' must match those of 'xvars' element of input 'object'.")}
        midx <- midx[-widx]
        newdata <- newdata[-widx]
        newDnames <- newDnames[-widx]
        midx <- match(newDnames,xnames)
      }
      lenmidx <- length(midx)
      if(lenmidx<nxvar){
        amidx <- 1:nxvar
        for(k in 1:lenmidx){amidx <- amidx[-which(amidx==midx[k])]}
        newDnames <- c(newDnames,xnames[amidx])
        newdata <- c(newdata,rep(NA,nxvar-lenmidx))
      }
      midx <- match(xnames,newDnames)
      newdata <- newdata[midx]
      names(newdata) <- xnames
    } # end if(is.null(newdata))
    
    ### check include
    newnames <- include
    oldnames <- object$tnames
    lnt <- length(newnames)
    midx <- match(newnames,oldnames)
    if(any(is.na(midx))){stop("Term names in 'include' must match those in 'tnames' element of input 'object'.")}
    newform <- as.formula(paste("y~",paste(newnames,collapse="+"),sep=""))
    trmfrm <- terms.formula(newform)
    et <- attr(trmfrm,"factors")
    etdim <- dim(et)
    etnames <- rownames(et)[2:etdim[1]]
    etidx <- match(etnames,xnames)
    Etab <- matrix(0L,nxvar,etdim[2])
    rownames(Etab) <- xnames
    Etab[etidx,] <- et[-1,]
    Etab <- rbind(0L,Etab)
    
    ### check for missing and transform
    newdim <- rep(NA,nxvar)
    for(k in 1:nxvar){
      if(sum(Etab[k+1L,])>0){
        if(any(is.na(newdata[[k]]))){
          stop(paste("Missing data in 'newdata' for",xnames[k],"term in 'include' input."))
        }
        if(any(object$type[[k]]==c("cub","cub0","per","tps"))){
          newdata[[k]] <- as.matrix(newdata[[k]]+0.0)
        } else {
          fidx <- match(newdata[[k]],object$modelspec$flvls[[k]])
          if(any(is.na(fidx))){stop(paste("Inappropriate 'newdata' for",xnames[k],"(factor levels don't match)."))}
          newdata[[k]] <- as.matrix(as.integer(fidx))
        }
        if(any(object$type[[k]]==c("cub","cub0","per"))){
          newdata[[k]] <- (newdata[[k]]-object$modelspec$xrng[[k]][1])/(object$modelspec$xrng[[k]][2]-object$modelspec$xrng[[k]][1])
        } 
        newdim[k] <- nrow(newdata[[k]])
        if(k>1L && newdim[k]!=max(newdim,na.rm=TRUE)){stop("Must have same number of observations for each covariate in newdata.")}
      } else {
        newdata[[k]]=NA
      } # end if(sum(Etab[k+1L,])>0)
    } # end for(k in 1:nxvar)
    
    ### make marginal reproducing kernel matrices
    rks <- makerkm(newdata,object$type,object$modelspec$myknots,object$modelspec$xrng,
                   pred=TRUE,tpsinfo=object$modelspec$tpsinfo)
    
    ### make design and penalty matrices
    if(lnt==length(oldnames) && effect=="all"){intid <- TRUE} else {intid <- includeint}
    dps <- ssadpm(newdata,object$type,rks[1:3],object$modelspec$myknots,Etab,pred=TRUE,effect,intid)
    Knames <- colnames(dps$Kmat)
    nknots <- nrow(object$modelspec$myknots[[1]])
    if(is.null(dps$Jmats)==FALSE){
      jdim <- dim(dps$Jmats)
      Jnames <- colnames(dps$Jmat)[seq(1,jdim[2],by=nknots)]
      Jlen <- length(Jnames)
    } else{Jlen <- NULL}
    
    ### get appropriate coefficients
    kidx <- jidx <- Jmat <- NULL
    pse <- NA
    gnames <- names(object$modelspec$gammas)
    nbf <- length(object$modelspec$coef)-nknots
    if(is.null(Knames)){
      if(is.null(dps$Jmats)){stop("No non/linear effect for terms in 'include' input.")}
      gamvec <- NULL
      for(j in 1:Jlen){
        xi <- strsplit(Jnames[j],":")
        xidx <- match(xi[[1]],gnames)
        gamvec <- c(gamvec,prod(object$modelspec$gammas[xidx]))
      }
      jidx <- (nbf+1):(nbf+nknots)
      Jmat <- dps$Jmats%*%kronecker(gamvec,diag(nknots))
      yhat <- Jmat%*%object$modelspec$coef[jidx]
      if(se.lp){pse <- sqrt(postvar(NULL,Jmat,object$modelspec$coef.csqrt[jidx,]))}
    } else if(is.null(dps$Jmats)){
      for(k in 1:nbf){kidx <- c(kidx,which(object$modelspec$Knames==Knames[k]))}
      kidx <- unique(kidx)
      yhat <- dps$Kmat%*%object$modelspec$coef[kidx]
      if(se.lp){pse <- sqrt(postvar(dps$Kmat,NULL,object$modelspec$coef.csqrt[kidx,]))}
    } else {
      for(k in 1:nbf){kidx <- c(kidx,which(object$modelspec$Knames==Knames[k]))}
      kidx <- unique(kidx)
      gamvec <- NULL
      for(j in 1:Jlen){
        xi <- strsplit(Jnames[j],":")
        xidx <- match(xi[[1]],gnames)
        gamvec <- c(gamvec,prod(object$modelspec$gammas[xidx]))
      }
      jidx <- (nbf+1):(nbf+nknots)
      Jmat <- dps$Jmats%*%kronecker(gamvec,diag(nknots))
      yhat <- cbind(dps$Kmat,Jmat)%*%object$modelspec$coef[c(kidx,jidx)]
      if(se.lp){pse <- sqrt(postvar(dps$Kmat,Jmat,object$modelspec$coef.csqrt[c(kidx,jidx),]))}
    }# end if(is.null(Knames))
    
    ### make fitted values
    if(object$family=="binomial"){
      fitvals <- 1/(1+exp(-yhat))
    } else if(object$family=="poisson"){
      fitvals <- exp(yhat)
    } else if(object$family=="Gamma"){
      yhat[yhat<=0] <- 10^-4
      fitvals <- 1/yhat
    } else if(object$family=="inverse.gaussian"){
      yhat[yhat<=0] <- 10^-4
      fitvals <- sqrt(1/yhat)
    } else if(object$family=="negbin"){
      fitvals <- exp(yhat)
    }
    
    ### collect new yhat
    pred <- list(fitted.values=as.numeric(fitvals),linear.predictors=as.numeric(yhat))
    if(se.lp){pred <- c(pred,list(se.lp=pse))} 
    if(design){pred <- c(pred,list(X=cbind(dps$Kmat,Jmat),ix=c(kidx,jidx)))}
    if(smoothMatrix){pred <- c(pred,list(S=tcrossprod(cbind(dps$Kmat,Jmat)%*%object$modelspec$coef.csqrt[c(kidx,jidx),])))}
    return(pred)
    
  }