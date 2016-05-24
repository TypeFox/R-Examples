summary.bce <- function(object,           # a bce-object, output of the function bce1() or BCE()
                        confInt=2/3, # confidence interval of values of composion matrix and ratio matrix
                        ...)            # additional arguments affecting the summary produced
  ## extract best, mean, sd, upper and lower boundaries, and covariance
  {
    if (!is.null(attributes(object)$A_not_null))
      {
        covariance <- cov(object$pars)
        
        w <- attributes(object)$A_not_null
        npig <- nrow(w)
        nalg <- ncol(w)
        outputlength <- nrow(object$pars)
        pignames <- attributes(object)$pignames
        algnames <- attributes(object)$algnames
        lw <- length(which(w))
        lp <- ncol(object$pars)
        Xratios <- attributes(object)$Xratios
        if (Xratios) nst <- (lp-lw)/(nalg-1) else nst <- (lp-lw)/nalg
        stnames <- attributes(object)$stnames
        
        mcmc.A <- array(0,dim=c(npig,nalg,outputlength),dimnames=list(pignames,algnames,NULL))
        mcmc.A[w] <- t(object$pars[,1:lw])
        mcmc.X <- array(dim=c(nalg,nst,outputlength),dimnames=list(algnames,stnames,NULL))
        if (Xratios)
          {
            suppressWarnings(Z <- matrix(c(1,-1,rep(0,nalg-1)),nalg,nalg-1))
            P <- c(rep(0,nalg-1),1)
            mcmc.Q <- array(dim=c(nalg-1,nst,outputlength),dimnames=list(NULL,stnames,NULL))
            mcmc.Q[] <- t(object$pars[,-(1:lw)])
            for (i in 1:outputlength) mcmc.X[,,i] <- Z%*%mcmc.Q[,,i]+P
          } else {
            mcmc.X[] <- t(object$pars[,-(1:lw)])
          }
        
        mean.pars <- apply(object$pars,2,mean)
        best.pars <- object$pars[which.min(object$SS),]
        sd.pars <- apply(object$pars,2,sd)
        last.pars <- object$pars[nrow(object$pars),]
        median.pars <- apply(object$pars,2,quantile,probs=1/2)
        ub.pars <- apply(object$pars,2,quantile,probs=(1+confInt)/2)
        lb.pars <- apply(object$pars,2,quantile,probs=(1-confInt)/2)
                         
        mean.A <- w; mean.A[w] <- mean.pars[1:lw]
        best.A <- w; best.A[w] <- best.pars[1:lw]
        sd.A <- w; sd.A[w] <- sd.pars[1:lw]
        last.A <- w; last.A[w] <- last.pars[1:lw]
        median.A <- w; median.A[w] <- median.pars[1:lw]
        ub.A <- w; ub.A[w] <- ub.pars[1:lw]
        lb.A <- w; lb.A[w] <- lb.pars[1:lw]

        if (Xratios)
          {
            suppressWarnings(Z <- matrix(c(1,-1,rep(0,nalg-1)),nalg,nalg-1))
            P <- c(rep(0,nalg-1),1)
            mean.Q <- matrix(mean.pars[-(1:lw)],nrow=nalg-1)
            mean.X <- Z%*%mean.Q+P
            best.Q <- matrix(best.pars[-(1:lw)],nrow=nalg-1)
            best.X <- Z%*%best.Q+P
            sd.X <- apply(mcmc.X,1:2,sd)
            last.Q <- matrix(last.pars[-(1:lw)],nrow=nalg-1)
            last.X <- Z%*%last.Q+P
            median.Q <- matrix(median.pars[-(1:lw)],nrow=nalg-1)
            median.X <- Z%*%median.Q+P
            ub.Q <- matrix(ub.pars[-(1:lw)],nrow=nalg-1)
            ub.X <- Z%*%ub.Q+P
            lb.Q <- matrix(lb.pars[-(1:lw)],nrow=nalg-1)
            lb.X <- Z%*%lb.Q+P
          } else {
            mean.X <- matrix(mean.pars[-(1:lw)],nrow=nalg)
            best.X <- matrix(mean.pars[-(1:lw)],nrow=nalg)
            sd.X <- matrix(sd.pars[-(1:lw)],nrow=nalg)
            last.X <- matrix(last.pars[-(1:lw)],nrow=nalg)
            median.X <- matrix(median.pars[-(1:lw)],nrow=nalg)
            ub.X <- matrix(ub.pars[-(1:lw)],nrow=nalg)
            lb.X <- matrix(lb.pars[-(1:lw)],nrow=nalg)
          }

        rownames(mean.A) <- rownames(best.A) <- rownames(sd.A) <- rownames(last.A) <- rownames(median.A) <- rownames(ub.A) <- rownames(lb.A) <- pignames
        colnames(mean.A) <- colnames(best.A) <- colnames(sd.A) <- colnames(last.A) <- colnames(median.A) <- colnames(ub.A) <- colnames(lb.A) <- 
          rownames(mean.X) <- rownames(best.X) <- rownames(sd.X) <- rownames(last.X) <- rownames(median.X) <- rownames(ub.X) <- rownames(lb.X) <- algnames
        colnames(mean.X) <- colnames(best.X) <- colnames(sd.X) <- colnames(last.X) <- colnames(median.X) <- colnames(ub.X) <- colnames(lb.X) <- stnames

        return(list(meanA=mean.A,
                    meanX=mean.X,
                    bestA=best.A,
                    bestX=best.X,
                    sdA=sd.A,
                    sdX=sd.X,
                    lastA=last.A,
                    lastX=last.X,
                    medianA=median.A,
                    medianX=median.X,
                    ubA=ub.A,
                    ubX=ub.X,
                    lbA=lb.A,
                    lbX=lb.X,
                    covar=covariance))
        
        
      } else {
        
        with(object,{

          nalg <- dim(Rat)[1]
          lr <- length(Rat)/length(logp)
          lx <- length(X)/length(logp)

          
          w <- which.max(logp)
          bestLogp <- logp[w]

          bestRat <- Rat[,,w]
          meanrat <- rowMeans(Rat,dims=2)

          quantile1 <- function(x) quantile(x,probs=c((1-confInt)/2,1/2,(1+confInt)/2))
          quantilerat <- apply(Rat,1:2,quantile1)
          lbrat <- quantilerat[1,,]
          ubrat <- quantilerat[3,,]
          sdrat <- apply(Rat,1:2,sd)

          if (is.matrix(X)) {
            firstX <- X[,1]
            bestX <- X[,w]
            meanX <- rowMeans(X)
            quantileX <- apply(X,1,quantile1)
            lbX <- quantileX[1,]
            ubX <- quantileX[3,]
            sdX <- apply(X,1,sd)
          } else{
            firstX <- X[,,1]
            bestX <- X[,,w]    
            meanX <- rowMeans(X,dims=2)
            quantileX <- apply(X,1:2,quantile1)
            lbX <- quantileX[1,,]
            ubX <- quantileX[3,,]
            sdX <- apply(X,1:2,sd)
          }

          bestDat <- bestX%*%bestRat

          if (all(sdrat==0)) covrat <- 0 else
          {
            covratnames <- vector(length=lr)
            for (i in 1:lr) covratnames[i] <- paste("Rat(",(i-1)%%nalg+1,",",(i-1)%/%nalg+1,")",sep="")
            covrat <- var(matrix(aperm(Rat,c(3,1,2)),ncol=lr,dimnames=list(NULL,covratnames))[,sdrat>1e-8],na.rm=TRUE)
          }

          covXnames <- vector(length=lx)
          for (i in 1:lx) covXnames[i] <- paste("x(",(i-1)%/%nalg+1,",",(i-1)%%nalg+1,")",sep="")
          covX <- var(matrix(aperm(X),ncol=lx,dimnames=list(NULL,covXnames)),na.rm=TRUE)
          
          ## output
          return(invisible(list(firstX=firstX,        # X determined through least squares regression from the initial ratio matrix and the data matrix
                                bestRat=bestRat,      # ratio matrix for which the posterior probability is maximal
                                bestX=bestX,          # composition matrix for which the posterior probability is maximal
                                bestLogp=bestLogp,    # maximal posterior probability
                                bestDat=bestDat,      # product of bestRat and bestX
                                meanRat=meanrat,      # means of the elements of the ratio matrix
                                sdRat=sdrat,          # standard deviation of the elements of the ratio matrix
                                lbRat=lbrat,          # lower boundary of the confidence interval of the elements of the ratio matrix
                                ubRat=ubrat,          # upper boundary of the confidence interval of the elements of the ratio matrix
                                covRat=covrat,        # covariance matrix of the elements of the ratio matrix
                                meanX=meanX,          # means of the elements of the composition matrix
                                sdX=sdX,              # standard deviation of the elements of the composition matrix
                                lbX=lbX,              # lower boundary of the confidence interval of the elements of the composition matrix
                                ubX=ubX,              # upper boundary of the confidence interval of the elements of the composition matrix
                                covX=covX             # covariance matrix of the elements of the composition matrix
                                )))
        })
      }
  } # end function summary.bce


