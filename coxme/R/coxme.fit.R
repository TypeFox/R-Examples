# Automatically generated from all.nw using noweb
coxme.fit <- function(x, y, strata, offset, ifixed, control,
                        weights, ties, rownames, 
                        imap, zmat, varlist, vparm, itheta,
                        ntheta, ncoef, refine.n, is.variance) {
    #     time0 <- proc.time() #debugging line
    n <-  nrow(y)
    if (length(x) ==0) nvar <-0
    else nvar <- ncol(as.matrix(x))
    
    if (missing(offset) || is.null(offset)) offset <- rep(0.0, n)
    if (missing(weights)|| is.null(weights))weights<- rep(1.0, n)
    else {
        if (any(weights<=0)) stop("Invalid weights, must be >0")
        }
    if (ncol(y) ==3) {
        if (length(strata) ==0) {
            sorted <- cbind(order(-y[,2], y[,3]), 
                            order(-y[,1]))
            newstrat <- n
            }
        else {
            sorted <- cbind(order(strata, -y[,2], y[,3]),
                            order(strata, -y[,1]))
            newstrat  <- cumsum(table(strata))
            }
        status <- y[,3]
        ofile <-  'agfit6b'
        rfile <-  'agfit6d'
        coxfitfun<- agreg.fit
        }
    else {
        if (length(strata) ==0) {
            sorted <- order(-y[,1], y[,2])
            newstrat <- n
            }
        else {
            sorted <- order(strata, -y[,1], y[,2])
            newstrat <-  cumsum(table(strata))
            }
        status <- y[,2]
        ofile <- 'coxfit6b' # fitting routine
        rfile <- 'coxfit6d' # refine.n routine
        coxfitfun <- coxph.fit
        }
    if (is.null(ifixed) ) {
        ifixed <- rep(0., ncol(x))
        if (length(ifixed) ==0) ifixed <- NULL  #agreg.fit didn't like numeric(0)
    }
    else if (length(ifixed) != ncol(x))
        stop("Wrong length for initial parameters of the fixed effects")

    if (length(itheta)==0) itemp <- 0 else itemp <- control$iter.max
    fit0 <- coxfitfun(x,y, strata=strata, 
                      offset=offset, init=ifixed, weights=weights,
                      method=ties, rownames=1:nrow(y),
                      control=coxph.control(iter.max=itemp))
    loglik0 <- fit0$loglik[length(fit0$loglik)]  # in case of no covariates  
    control$inner.iter <- eval(control$inner.iter)
    kfun <- function(theta, varlist, vparm, ntheta, ncoef) {
        nrandom <- length(varlist)
        sindex <- rep(1:nrandom, ntheta) #which thetas to which terms

        tmat <- varlist[[1]]$generate(theta[sindex==1], vparm[[1]]) 
        dd <- dim(tmat)
        if (length(dd) !=2 || any(dd != rep(ncoef[1,1]+ncoef[1,2], 2)))
            stop("Incorrect dimensions for generated penalty matrix, term 1")
        if (!inherits(tmat, 'bdsmatrix')) 
            tmat <- bdsmatrix(blocksize=integer(0), blocks=numeric(0), rmat=tmat)
        if (nrandom ==1) return(tmat)

        # Need to build up the matrix by pasting up a composite R
        nsparse <- sum(tmat@blocksize)
        nrow.R <- sum(ncoef)
        ncol.R <- nrow.R - nsparse
        R <- matrix(0., nrow.R, ncol.R)
        indx1 <- 0                  #current column offset wrt intercepts
        indx2 <- sum(ncoef[,1]) -nsparse #current col offset wrt filling in slopes
        
        if (ncol(tmat) > nsparse) { #first matrix has an rmat component
            if (ncoef[1,1] > nsparse) { #intercept contribution to rmat
                irow <- 1:ncoef[1,1]  #rows for intercepts
                j <- ncoef[1,1] - nsparse   #number of dense intercept columns
                R[irow, 1:j] <- tmat@rmat[irow,1:j]
                indx1 <- j  #number of intercept processed so far
                
                if (ncoef[1,2] >0) {
                    # T[1-62, 3-66] of the example above
                    k <- 1:ncoef[1,2]
                    R[irow, k+indx2-nsparse] <- tmat@rmat[irow, k+j]
                }
                }
            else j <- 0
            
            if (ncoef[1,2] >0) { #has a slope contribution to rmat
                # T[63-128, 3-66] of the example above
                k <- 1:ncoef[1,2]
                R[k+indx2 +nsparse, k+ indx2] <- tmat@rmat[k+indx1, j+k]
                indx2 <- indx2 + ncoef[1,2] #non intercetps so far
                }
            }
     
    for (i in 2:nrandom) {
            temp <- as.matrix(varlist[[i]]$generate(theta[sindex==i], vparm[[i]]))
            if (any(dim(temp) != rep(ncoef[i,1]+ncoef[i,2], 2)))
                stop(paste("Invalid dimension for generated penalty matrix, term",
                           i))
            
            if (ncoef[i,1] >0)  { # intercept contribution
                #U or V [1-8, 1-8] in the example above
                j <- ncoef[i,1]
                R[indx1 +1:j + nsparse, indx1 +1:j] <- temp[1:j,1:j]
                
                if (ncoef[i,2] >0) {
                    # V[1-8, 9-24] in the example
                    k <- 1:ncoef[i,2]
                    R[indx1+ 1:j + nsparse, indx2 +k] <- temp[1:j, k+ j]
                    # V[9-24, 9-24]
                    R[indx2+k +nsparse, indx2 +k] <- temp[k+j, k+j]
                    }
                }
            else if (ncoef[i,2]>0) {
                k <- 1:ncoef[i,2]
                R[indx2+k +nsparse, indx2+k] <- temp
                }
            indx1 <- indx1 + ncoef[i,1]
            indx2 <- indx2 + ncoef[i,2]
            }
        bdsmatrix(blocksize=tmat@blocksize, blocks=tmat@blocks, rmat=R)
        }    
    if (length(itheta) >0) theta <- sapply(itheta, function(x) x[1])
    else theta <- numeric(0)
    dummy <- kfun(theta, varlist, vparm, ntheta, ncoef)
    if (is.null(dummy@rmat)) rcol <- 0
        else                 rcol <- ncol(dummy@rmat)
    npenal <- ncol(dummy)  #total number of penalized terms

    if (ncol(imap)>0) {
        findex <- matrix(0, nrow=sum(ncoef), ncol=ncol(imap))
        for (i in 1:ncol(imap)) findex[cbind(imap[,i], i)] <- 1
        }
    else findex <- 0  # dummy value

    if (is.null(control$sparse.calc)) {
        nevent <- sum(y[,ncol(y)])
        if (length(dummy@blocksize)<=1) nsparse<- 0
        else nsparse <- sum(dummy@blocksize)
        itemp <- max(c(0,imap)) - nsparse  #number of non-sparse intercepts
        
        if ((2*n) > (nevent*(nsparse-itemp))) control$sparse.calc <- 0
        else control$sparse.calc <- 1
        }

    ifit <- .C('coxfit6a', 
                   as.integer(n),
                   as.integer(nvar),
                   as.integer(ncol(y)),
                   as.double(c(y)),
                   as.double(cbind(zmat,x)),
                   as.double(offset),
                   as.double(weights),
                   as.integer(length(newstrat)),
                   as.integer(newstrat),
                   as.integer(sorted-1),
                   as.integer(ncol(imap)),
                   as.integer(imap-1),
                   as.integer(findex),
                   as.integer(length(dummy@blocksize)),
                   as.integer(dummy@blocksize),
                   as.integer(rcol),
                   means = double(nvar),
                   scale = double(nvar),
                   as.integer(ties=='efron'),
                   as.double(control$toler.chol),
                   as.double(control$eps),
                   as.integer(control$sparse.calc))
    means   <- ifit$means
    scale   <- ifit$scale
    logfun <- function(theta, varlist, vparm, kfun, ntheta, ncoef, 
                       init, fit0, iter, ofile) {
        gkmat <- gchol(kfun(theta, varlist, vparm, ntheta, ncoef))
        if (is.variance) {
            ikmat <- solve(gkmat)  #inverse of kmat, which is the penalty
            if (any(diag(ikmat) <=0)) { #Not an spd matrix
                return(0)  # return a "worse than null" fit
            }
            fit <- .C(ofile,
                      iter= as.integer(c(iter,iter)),
                      beta = as.double(init),
                      loglik = double(2),
                      as.double(ikmat@blocks),
                      as.double(ikmat@rmat),
                      hdet = double(1))
            ilik <- fit$loglik[2] -
                .5*(sum(log(diag(gkmat))) + fit$hdet)
        }
        else {
            # The variance functions have returned the inverse matrix
            fit <- .C(ofile,
                      iter= as.integer(c(iter,iter)),
                      beta = as.double(init),
                      loglik = double(2),
                      as.double(gkmat@blocks),
                      as.double(gkmat@rmat),
                      hdet = double(1))
            ilik <- fit$loglik[2] +
                .5*(sum(log(diag(gkmat))) - fit$hdet)
        }
        -(1+ ilik - fit0)
        }

    ishrink <- 0.7  # arbitrary guess
    init.coef <- c(rep(0., npenal), scale*fit0$coef* ishrink)

    if (length(itheta)==0) iter <- c(0,0)
    else {
        nstart <- sapply(itheta, length)
        if (all(nstart==1)) theta <- unlist(itheta)  #one starting guess
        else {
            #make a matrix of all possible starting estimtes
            testvals <- do.call(expand.grid, itheta)
            bestlog <- NULL
            for (i in 1:nrow(testvals)) {
                ll <- logfun(as.numeric(testvals[i,]), 
                             varlist, vparm, kfun, ntheta, ncoef, 
                             init=init.coef, loglik0, 
                             control$inner.iter, ofile)
                if (is.finite(ll)) {
                    #ll calc can fail if someone picks a very bad starting guess
                    if (is.null(bestlog) || ll < bestlog) {  
                        # (optim is set up to minimize)
                        bestlog <- ll
                        theta <- as.numeric(testvals[i,])
                    }
                }
            }
            if (is.null(bestlog))
                stop("No starting estimate was successful")
        }
        # Finally do the fit
        logpar <- list(varlist=varlist, vparm=vparm, 
                       ntheta=ntheta, ncoef=ncoef, kfun=kfun,
                       init=init.coef, fit0= loglik0,
                       iter=control$inner.iter,
                       ofile=ofile)
        mfit <- do.call('optim', c(list(par= theta, fn=logfun, gr=NULL), 
                               control$optpar, logpar))
        theta <- mfit$par
        iter <- mfit$counts[1] * c(1, control$inner.iter)
    }
    gkmat <- gchol(kfun(theta, varlist, vparm, ntheta, ncoef))
    if (is.variance) {
        ikmat <- solve(gkmat)  #inverse of kmat, which is the penalty
        fit <- .C(ofile,
                  iter= as.integer(c(0, control$iter.max)),
                  beta = as.double(c(rep(0., npenal), fit0$coef*scale)),
                  loglik = double(2),
                  as.double(ikmat@blocks),
                  as.double(c(ikmat@rmat,0)),
                  hdet = double(1))
        ilik <- fit$loglik[2] -
            .5*(sum(log(diag(gkmat))) + fit$hdet)
    } else {
        fit <- .C(ofile,
                  iter= as.integer(c(0, control$iter.max)),
                  beta = as.double(c(rep(0., npenal), fit0$coef*scale)),
                  loglik = double(2),
                  as.double(gkmat@blocks),
                  as.double(c(gkmat@rmat,0)),
                  hdet = double(1))
        ilik <- fit$loglik[2] +
            .5*(sum(log(diag(gkmat))) - fit$hdet)
    }
    iter[2] <- iter[2] + fit$iter[2]
    nfrail <- nrow(ikmat)  #total number of penalized terms
    nsparse <- sum(ikmat@blocksize)
    nvar2  <- nvar + (nfrail - nsparse)  # total number of non-sparse coefs
    nvar3  <- as.integer(nvar + nfrail)  # total number of coefficients
    btot   <- length(ikmat@blocks)

    fit3 <- .C('coxfit6c',
                   u    = double(nvar3),
                   h.b  = double(btot),
                   h.r  = double(nvar2*nvar3),
                   hi.b = double(btot),
                   hi.r = double(nvar2*nvar3),
                   hrank= integer(1),
                   as.integer(ncol(y))
                   )
        if (nvar2 ==0) {
            hmat <- new('gchol.bdsmatrix', Dim=c(nvar3, nvar3),
                        blocksize=ikmat@blocksize, blocks=fit3$h.b,
                        rmat=matrix(0,0,0), rank=fit3$hrank,
                        Dimnames=list(NULL, NULL))
            hinv <- bdsmatrix(blocksize=ikmat@blocksize, blocks=fit3$hi.b)
            }
        else {
            rmat1 <- matrix(fit3$h.r, nrow=nvar3)
            rmat2 <- matrix(fit3$hi.r, nrow=nvar3)
            if (nvar ==1 ) {
                rmat1[nvar3,] <- rmat1[nvar3,]/scale
                rmat2[nvar3,] <- rmat2[nvar3,]/scale
                rmat1[,nvar2] <- rmat1[,nvar2]*scale
                rmat2[,nvar2] <- rmat2[,nvar2]/scale
                rmat1[nvar3,nvar2] <- rmat1[nvar3,nvar2]*scale^2
                u <- fit3$u  # the efficient score vector U
                u[nvar3] <- u[nvar3]*scale
                }
            else if (nvar >1) {
                temp <- seq(to=nvar3, length=length(scale))
                u <- fit3$u
                u[temp] <- u[temp]*scale
                rmat1[temp,] <- (1/scale)*rmat1[temp,] #multiply rows* scale
                rmat2[temp,] <- (1/scale)*rmat2[temp,] 

                temp <- temp-nsparse          #multiply cols
                rmat1[,temp] <- rmat1[,temp] %*% diag(scale)
                rmat2[,temp] <- rmat2[,temp] %*% diag(1/scale)
                temp <- seq(length=length(scale), to=length(rmat1), by=1+nvar3)
                rmat1[temp] <- rmat1[temp]*(scale^2)    #fix the diagonal
                }
            hmat <- new('gchol.bdsmatrix', Dim=c(nvar3, nvar3),
                        blocksize=ikmat@blocksize, blocks=fit3$h.b,
                        rmat= rmat1, rank=fit3$hrank,
                        Dimnames=list(NULL, NULL))
            hinv <- bdsmatrix(blocksize=ikmat@blocksize, blocks=fit3$hi.b,
                              rmat=rmat2)
            }
    traceprod <- function(H, P) {
        #block-diagonal portions will  match in shape
        nfrail <- nrow(P)  #penalty matrix
        nsparse <- sum(P@blocksize)
        if (nsparse >0) {
            temp1 <- 2*sum(H@blocks * P@blocks) -
                     sum(diag(H)[1:nsparse] * diag(P)[1:nsparse])
            }
        else temp1 <- 0
        
        if (length(P@rmat) >0) {
            #I only want the penalized part of H
            rd <- dim(P@rmat)
            temp1 <- temp1 + sum(H@rmat[1:rd[1], 1:rd[2]] * P@rmat)
            }
        temp1
        }

    df <- nvar + (nfrail - traceprod(hinv, ikmat))
    if (refine.n > 0) {
        rdf <- control$refine.df
        nfrail <- ncol(gkmat)
        hmatb <- hmat[1:nfrail, 1:nfrail]
        if (control$refine.method == "control") {
            #create the random t-variate with variance H-inverse
            bmat <- matrix(rnorm(nfrail*refine.n), ncol=refine.n)
            bmat <- backsolve(hmatb, bmat, upper=TRUE) /
                rep(sqrt(rchisq(refine.n, df=rdf)/rdf), each=nfrail)
            bmat2 <- bmat + fit$beta[1:nfrail]  #recenter
        }
        else if (control$refine.method == "direct") {
            bmat <- matrix(rnorm(nfrail*refine.n), ncol=refine.n)
            bmat2 <- gkmat %*% bmat
        }
        else stop("Unrecognized value for refine.method")
        
        rfit <- .C(rfile,
                   as.integer(refine.n),
                   as.double(fit$beta),
                   as.double(bmat2),
                   loglik = double(refine.n), dup=FALSE)

        if (control$refine.method == "direct") {
            temp <- max(rfit$loglik)   #keep exp() in range
            errhat <- exp(rfit$loglik - temp) 
            mtemp <- mean(errhat)             #estimated integral
            stemp <- sqrt(var(errhat)/refine.n)   #std of the estimate
            r.correct <- c(correction = log(mtemp) + temp -ilik, std= stemp/mtemp)
        }
        else {
            # Penalty terms
            penalty1 <- colSums(bmat2*(ikmat %*% bmat2))/2
            penalty2 <- rowSums((t(bmat) %*% hmatb)^2 )/2

            # Constant for the Gaussian density,  and density of the t-dist (logs)
            gdens <- -0.5* (sum(log(diag(gkmat))) + nfrail*log(2* pi))
            logdet <- -sum(log(diag(hmatb)))
            tdens <- lgamma((nfrail + rdf)/2) - 
                (lgamma(rdf/2) + 0.5*(logdet + nfrail* log(pi*rdf) + 
                                      (nfrail+ rdf)* log(1 + 2*penalty2/rdf)))
        
            # Add it up, we have to be very careful about round-off
            n1 <- rfit$loglik + gdens - (penalty1 + ilik + tdens)
            n2 <- fit$loglik[2] + gdens - (penalty2 + ilik + tdens)
            temp <- max(n1, n2)  #scale so the largest value is about 1
            errhat <- (exp(n1-temp) - exp(n2-temp)) * exp(temp)
            #errhat <- (exp(rfit$loglik -(penalty1 + ilik)) - 
            #        exp(fit$loglik[2]- (penalty2 + ilik))) * exp(gdens-tdens)
      
            mtemp <- mean(errhat)             #estimated integral
            stemp <- sqrt(var(errhat)/refine.n)   #std of the estimate
            r.correct <- c(correction= log(1+ mtemp), std=stemp/(1 +mtemp)) 
        }
    }
    .C("coxfit6e", as.integer(ncol(y)))  #release memory
    idf <- nvar + sum(ntheta)
    fcoef <- fit$beta[1:nfrail]
    penalty <- sum(fcoef * (ikmat %*% fcoef))/2

    if (nvar > 0) {
        out <- list(coefficients = fit$beta[-(1:nfrail)]/scale, frail=fcoef, 
             theta=theta, penalty=penalty,
             loglik=c(fit0$log[1], ilik, fit$log[2]), variance=hinv,
             df=c(idf, df), hmat=hmat, iter=iter, control=control,
             u=u, means=means, scale=scale)
        }
    else out <- list(coefficients=NULL, frail=fcoef, 
                     theta=theta, penalty=penalty,
              loglik=c(fit0$log[1], ilik, fit$log[2]), variance=hinv,
              df=c(idf, df), hmat=hmat, iter=iter, control=control,
              u=fit3$u, means=means, scale=scale)    

    if (refine.n>0) {
        out$refine <- r.correct
        #The next line can be turned on for detailed tests in refine.R
        #  The feature is not documented in the manual pages, only
        #  here.
        if (control$refine.detail) {
            if (control$refine.method== "control")
                out$refine.detail <-list(loglik=rfit$loglik, bmat=bmat2, 
                                         tdens=tdens,
                                         penalty1=penalty1, penalty2=penalty2,
                                         gdens=gdens, errhat=errhat, gkmat=gkmat)
            else out$refine.detail <- list(loglik=rfit$loglik, bmat=bmat2,
                                           errhat=errhat, gkmat=gkmat)
        }
    }
    out
}
