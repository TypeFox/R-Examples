# Automatically generated from all.nw using noweb
coxmeFull <- function(collapse=FALSE) {
    collapse <- collapse
    # Because of R's lexical scoping, the values of the options
    #  above, at the time the line below is run, are known to the
    #  initialize function
    initialize <- function(vinit, fixed, intercept, G, X,  control) {
        ngroup <- min(length(G), ncol(G))
        nvar   <- min(length(X), ncol(X))  # a NULL or a nx0 matrix yields 0
        sparse <- control$sparse
        initmatch <- function(namelist, init) {
            if (is.null(names(init))) iname <- rep('', length(init))
            else iname <- names(init)
            
            indx1 <- pmatch(iname, namelist, nomatch=0, duplicates.ok=TRUE)
            if (any(iname=='')) {
                temp <- 1:length(namelist)
                if (any(indx1>0)) temp <- temp[-indx1]   #already used
                indx2 <- which(iname=='')
                n <- min(length(indx2), length(temp))
                if (n>0) indx1[indx2[1:n]] <- temp[1:n]
                }
            indx1
            }
        
        if (ngroup==0) {
            if (intercept)
                return(list(error=("Invalid null random term (1|1)")))
            else {
                xname <- dimnames(X)[[2]]
                if (length(vinit) >0) {
                  temp <- initmatch(xname[1], vinit)
                  if (any(temp==0)) 
                      return(list(error=paste('Element', which(temp==0),
                                              'of initial values not matched')))
                  else {
                      if (vinit <=0) return(list(error="Invalid variance value, must be >0")) 
                      theta <- vinit
                      }
                  }
                else theta <- control$varinit *.5 / mean(sqrt(apply(X,2,var)))
                  
                if (length(fixed) >0) {
                    temp <- initmatch(xname[1], fixed)
                    if (any(temp==0))
                        return(list(error=paste('Element', which(temp==0),
                                                'of fixed variance values not matched')))
                    else theta <- fixed
                    which.fixed <- TRUE
                    if (theta <=0) return(list(error="Invalid variance value, must be >0"))
                    }
                else which.fixed <- FALSE

                xmap <- matrix(0L, nrow=nrow(X), ncol=ncol(X))
                for (i in 1:ncol(X)) xmap[,i] <- i

                list(theta=list(log(theta))[!which.fixed], imap=NULL, X=X, xmap=xmap,
                         parms=list(fixed=which.fixed, theta=theta[1],
                                    xname=xname, case=1))
                }
            }
        else {
            if (nvar==0) {
                gname <-  names(G)
                ntheta <- length(gname)
                itheta <- vector('list', length=ntheta)
                for (i in 1:ntheta) itheta[[i]] <- control$varinit
                if (ntheta >1) {
                    for (i in 2:ntheta) gname[i] <- paste(gname[i-1], gname[i], sep='/')
                    gname <- rev(gname) 
                    }
                names(itheta) <- gname

                if (length(vinit) >0) {
                    temp <- initmatch(gname, vinit)
                    if (any(temp==0))
                        return(list(error=paste('Element', which(temp==0),
                                                'of initial values not matched')))
                    else itheta[temp] <- vinit
                    if (any(unlist(vinit) <=0))
                        return(list(error='Invalid initial value'))
                    }

                theta <- rep(0, ntheta)   # the filler value does not matter
                which.fixed <- rep(FALSE, ntheta)
                if (length(fixed)>0) {
                    temp <- initmatch(gname, fixed)
                    if (any(temp==0))
                        return(list(error=paste('Element', which(temp==0),
                                                 'of variance values not matched')))
                    else theta[temp] <- unlist(fixed)
                    which.fixed[temp] <- TRUE
                    }
                }

            # Deal with random slope terms
            if (ngroup ==1) {
                gtemp <- as.factor(G[[1]])[,drop=TRUE] #drop unused levels
                nlevel <- length(levels(gtemp))
                gfrac <- table(gtemp)/ length(gtemp)
                if (nlevel >= sparse[1] && any(gfrac <= sparse[2])) {
                    indx <- order((gfrac> sparse[2]), 1:nlevel)  #False then True for order
                    nsparse <- sum(gfrac <= sparse[2])
                    if (nsparse== nlevel) vmat<- bdsI(nsparse)
                    else {
                        k <- nlevel - nsparse  #number of non-sparse levels
                        rmat <- matrix(0., nrow=nlevel, ncol=k)
                        rmat[seq(nsparse+1, by= nlevel+1, length=k)] <- 1.0
                        vmat <- bdsmatrix(blocksize=rep(1,nsparse), 
                                          blocks= rep(1,nsparse), rmat=rmat)
                        }
                    }
                else {
                    vmat <- diag(nlevel)
                    indx <- 1:nlevel
                    nsparse <- 0
                    }
                imap <- matrix(match(as.numeric(gtemp), indx))
                levellist <- list((levels(gtemp))[indx]) 
                varlist <- list(vmat)
                if (nvar==0) 
                    return(list(imap=imap, X=NULL, 
                                theta=(lapply(itheta, log))[!which.fixed], 
                                parms=list(varlist=varlist, theta=theta, levellist=levellist,
                                           fixed=which.fixed, case=2, gname=gname)))
                    }
            else {
                if (collapse) {
                    gtemp <- expand.nested(G)
                    n <- nrow(G)

                    #Sparse?
                    gfrac <- table(gtemp[,ngroup])/ n
                    nlevel <- length(levels(gtemp))
                    if (nlevel > sparse[1] && any(gfrac <= sparse[2])) {
                            is.sparse <- (gfrac <= sparse[2])[as.numeric(gtemp[,ngroup])]
                            block.sparse <- tapply(is.sparse, G[,1], all)
                            nsparse
                            }
                    else block.sparse <- 0

                    if (sum(block.sparse > 1)) { #sparse blocks exist, make them list first
                        border <- order(!block.sparse, 1:nlevel)
                        G[,1] <- factor(gtemp[,1], levels=levels(G[,1])[border])
                        G <- expand.nested(G)
                        }

                    G <- rev(expand.nested(G))
                    levellist <- lapply(G, levels)
                    nlevel <-  sapply(levellist, length)
                    imap <- matrix(as.numeric(G[,1]))

                    varlist <- vector('list', ngroup)
                    indx <- match(levellist[[1]], G[[1]])  #an ordered set of unique rows
                    for (i in 1:ngroup) 
                        varlist[[i]] <- bdsBlock(1:nlevel[1], G[indx,i])

                    if (sum(block.sparse) <=1) {#make them all ordinary matrices
                        for (i in 1:ngroup) varlist[[i]] <- as.matrix(varlist[[i]])
                        }
                    else {
                        if (!all(block.sparse)) { # Only a part is sparse
                            tsize <- sum(temp@blocksize[1:sum(block.sparse)]) # sparse coefs
                            sparse <- 1:tsize  #sparse portion, remainder is dense
                            smat <- (varlist[[ngroup]])[1:sparse, 1:sparse]
                            varlist[[ngroup]]
                            rmat <- matrix(0, sum(tsize), nlevel[1])
                            rmat[seq(by=nrow(rmat)+1, to=length(rmat), length=ncol(rmat))] <- 1.0
                            varlist[[ngroup]] <- bdsmatrix(blocksize=smat@blocksize,
                                                           blocks=smat@block, rmat=rmat)
                            } 
                        varlist <- bdsmatrix.reconcile(varlist)
                        }

                    if (nvar==0) {
                        return(list(imap=matrix(as.numeric(G[,1])), X=NULL, 
                                    theta=lapply(itheta, log)[!which.fixed],
                                    parms=list(varlist=varlist, theta=theta, 
                                               fixed=which.fixed, gname=names(G),
                                               levellist=levellist, case=3, collapse=TRUE)))
                        }
                    }
                else {
                    G <- rev(expand.nested(G))  #the last shall be first
                    imap <- matrix(0L, nrow=nrow(G), ncol=ngroup)
                    imap[,1] <- as.numeric(G[,1])
                    for (i in 2:ngroup) 
                        imap[,i] <- as.numeric(G[,i]) + max(imap[,i-1])
                    levellist <- lapply(G, levels)
                    nlevel <- sapply(levellist, length)

                    # Sparsity?
                    gtemp <- G[,1]
                    gfrac <- table(gtemp)/ length(gtemp)
                    if (nlevel[1] > sparse[1] && any(gfrac <= sparse[2])) {
                        indx <- order((gfrac> sparse[2]), 1:nlevel[1])
                        nsparse <- sum(gfrac <= sparse[2])

                        imap[,1] <- match(as.integer(gtemp), indx) 
                        levellist[[1]] <- (levellist[[1]])[indx]
                        }
                    else  nsparse <- 0  #a single sparse element is the same as dense
                    if (nsparse==0) tmat <- diag(sum(nlevel))
                    else tmat <- bdsmatrix(blocksize=rep(1L, nsparse), blocks=rep(1., nsparse),
                                       rmat=matrix(0., nrow=sum(nlevel), ncol=sum(nlevel)-nsparse))
                    varlist <- vector('list', ngroup) 
                    for (i in 1:ngroup) {
                        temp <- rep(0., nrow(tmat))
                        temp[unique(imap[,i])] <- 1.0
                        temp2 <- tmat
                        diag(temp2) <- temp
                        varlist[[i]] <- temp2
                        }

                    if (nvar==0) 
                        return(list(imap=imap, X=NULL, 
                                    theta=(lapply(itheta, log))[!which.fixed],
                                    parms=list(nlevel=nlevel, varlist=varlist, gname=names(G),
                                               fixed=which.fixed, levellist=levellist, 
                                               theta=theta, case=3, collapse=FALSE)))           
                    }
                }

            #Deal with slopes
            if (nvar > 0) {
                xvar  <- apply(X,2,var)
                cordefault <- control$corinit
                itheta <- list()
                is.variance <- NULL
                if (intercept) {
                    itheta <- c(itheta, list(control$varinit))
                    for (i in 1:ncol(X)) itheta <- c(itheta, list(control$corinit)) #correlations
                    is.variance <- c(TRUE, rep(FALSE, ncol(X)))
                    }
                for (i in 1:ncol(X)) {
                    itheta <- c(itheta, list(control$varinit/xvar))  #variance
                    if (i < ncol(X)) {
                        for (j in (i+1):ncol(X))  itheta <- c(itheta, list(cordefault))
                        is.variance <- c(is.variance, TRUE, rep(FALSE, ncol(X)-i))
                        }
                    else is.variance <- c(is.variance, TRUE)
                    }
                itheta <- rep(itheta, ngroup)
                is.variance <- rep(is.variance, ngroup)    
                    
                xname <- dimnames(X)[[2]]
                name.temp <- outer(xname, xname, paste, sep=":")
                diag(name.temp) <- xname
                name.temp <- name.temp[row(name.temp) >= col(name.temp)]
                tname <- ""
                gname <- names(G)
                for (i in 1:ngroup) {
                    if (intercept)
                        tname <- c(tname, gname[i], paste(xname, gname[i], sep=':'),
                                   paste(name.temp, gname[i], sep='/'))
                    else tname <- paste(name.temp, gname[i], sep='/')
                    }

                # Now replace selected values with the user's input
                if (length(vinit) > 0) {
                    temp <- initmatch(tname, vinit)
                    if (any(temp==0))
                        return(list(error=paste('Element(s)', which(temp==0),
                                                'of initial values not matched')))
                    else itheta[temp] <- unlist(vinit)
                    }
                              
                which.fixed <- rep(FALSE, length(itheta))
                if (length(fixed) > 0) {
                    temp <- initmatch(tname, fixed)
                    if (any(temp==0))
                      return(list(error=paste('Element(s)', which(temp==0),
                                              'of fixed variance values not matched')))
                    else itheta[temp] <- unlist(fixed)
                    which.fixed[temp] <- TRUE
                    }

                # Check for legality of the values
                tmat <- diag(nvar+ intercept)  #dummy variance/cov matrix
                tmat <- tmat[row(tmat) >= col(tmat)]
                vindx <- which(tmat==1)  #indices of the variance terms within each group
                cindx <- which(tmat==0)  #indices of the correlations

                if (any(unlist(itheta[is.variance]) <=0)) 
                        return(list(error="Variances must be >0"))
                if (any(unlist(itheta[!is.variance]) <=-1) || any(unlist(itheta[!is.variance]) >=1))
                        return(list(error="Correlations must be between 0 and 1"))
                xnew <- matrix(0., nrow=nrow(X), ncol=nvar*ncol(imap))
                xmap <- matrix(0L, nrow=nrow(X), ncol=ncol(X)*ncol(imap))
                xoffset <- (intercept)* max(imap)
                k <- 1
                for (j in 1:nvar) {
                    for (i in 1:ncol(imap)) { 
                        xnew[,k] <- X[,j]
                        xmap[,k] <- imap[,i] + xoffset
                        k <- k+1
                        xoffset <- xoffset + max(imap)
                        }
                    }

                # Transform correlations to (1+rho)/(1-rho) scale, which is used for the saved
                #  parameters, and make a copy.  The initial parameters are then log transformed
                itheta[!is.variance] <- lapply(itheta[!is.variance], function(rho) (1+rho)/(1-rho))
                theta <- sapply(itheta, function(x) x[1])

                itheta <- lapply(itheta, log)

                if (intercept)
                    list(theta=itheta[!which.fixed], imap=imap, X=xnew, xmap=xmap, 
                             parms=list(theta=theta, fixed=which.fixed, 
                                        nlevel=nlevel, levellist=levellist,
                                        nvar=nvar, gname=names(G), varlist=varlist,
                                        xname=dimnames(X)[[2]], intercept=intercept,
                                        xname=dimnames(X)[[2]], case=4, collapse=collapse))
                else list(theta=itheta[!which.fixed], imap=NULL, X=xnew, xmap=xmap, 
                             parms=list(theta=theta, fixed=which.fixed, 
                                        nlevel=nlevel, levellist=levellist,
                                        nvar=nvar, gname=names(G), varlist=varlist,
                                        xname=dimnames(X)[[2]], intercept=intercept,
                                        xname=dimnames(X)[[2]], case=4, collapse=collapse))
                }
            }
        }
    generate= function(newtheta, parms) {
        theta <- parms$theta
        if (length(newtheta)>0) theta[!parms$fixed] <- 
            exp(pmax(-36, pmin(36, newtheta)))

        if (parms$case==1) return(diag(length(parms$xname)) * theta)

        if (parms$case==2) return(theta*parms$varlist[[1]])
        if (parms$case==3) {
            temp <- theta[1]* parms$varlist[[1]]
            for (i in 2:length(parms$varlist)) 
                temp <- temp + theta[i]*parms$varlist[[i]]
            return(temp)
            }
        if (parms$case==4) {
            ngroup <- length(parms$nlevel)
            n1 <- sum(parms$nlevel)          #number of intercept coefs
            nvar <- parms$nvar               #number of covariates
            n2 <-   n1*nvar                  #number of slope coefs
            theta.per.group <- length(theta)/ngroup
            tindx <- seq(1, by=theta.per.group, length=ngroup)
            
            addup <- function(theta, p=parms) {
                tmat <- theta[1] * p$varlist[[1]]
                if (length(theta) >1) {
                    for (i in 2:length(theta)) tmat <- tmat + theta[i]*p$varlist[[i]]
                    }
                tmat
                }
            
            if (parms$intercept) {
                # upper left corner (has no covarinaces)
                ivar <- theta[tindx]  #variances of the intercepts
                corner <- addup(ivar)
                if (inherits(corner, 'bdsmatrix')) {
                    nsparse <- sum(corner@blocksize)
                    rmat <- matrix(0., nrow=n1+n2, ncol=n1+n2 - nsparse)
                    if (nsparse < n1) rmat[1:n1, 1:(n1-nsparse)] <- corner@rmat
                    }
                else {
                    nsparse <- 0
                    rmat <- matrix(0., n1+n2, n1+n2)
                    rmat[1:n1, 1:n1] <- corner
                    }

                # Covariances with the intercept
                for (i in 1:nvar) {
                    xvar <- theta[i+nvar+tindx]  #variance of the slope
                    xcor <- (theta[i+tindx]-1)/(theta[i+tindx]+1)  # correlation
                    icov <- xcor * sqrt(xvar * ivar)               # covariance
                    rmat[1:n1, 1:n1 +i*n1 -nsparse] <- as.matrix(addup(icov)) 
                    }
                irow <- n1
                icol <- n1 - nsparse
                theta <- theta[-(1:(1+nvar))]  #these thetas are 'used up'
                }
            else {
                irow <- icol <- 0
                rmat <- matrix(0., n2,n2)
                }

            # covariates
            offset1 <- 0
            for (i in 1:nvar) {
                xvar <- theta[offset1 + tindx]  #variance of the slope
                rmat[1:n1 + irow, 1:n1 + icol] <- as.matrix(addup(xvar))
                # covariate-covariate
                if (i<nvar) {
                    offset2 <- offset1 + 1 + nvar-1
                    for (j in (i+1):nvar) {
                        icol <- irow + n1
                        zvar <- theta[offset2 + tindx]-1
                        zcor <- (theta[j+offset2+tindx] -1)/(theta[j+offset2 +tindx]+1)
                        zcov <- sqrt(xvar*zvar) * zcor
                        rmat[1:n1+ irow, 1:n1 + icol] <- as.matrix(addup(zcov))
                        offset2 <- offset2 + 1 + nvar -j
                        }
                    offset1 <- offset1 + 1 + nvar- i
                    icol <- irow <- irow+n1
                    }
                }

            if (parms$intercept && inherits(corner, 'bdsmatrix'))
                bdsmatrix(blocksize=corner@blocksize, blocks=corner@blocks, rmat=rmat)
            else bdsmatrix(blocksize=integer(0), blocks=numeric(0), rmat=rmat)
            }
        }
    wrapup <- function(theta, b, parms) {
        newtheta <- parms$theta
        if (length(theta)) newtheta[!parms$fixed] <- exp(theta)
        
        if (parms$case==1) {
            theta <- list(c('(Shrinkage)' = newtheta[1]))
            names(theta) <- '1'
            names(b) <- parms$xname
            return(list(theta =theta, b=list('1'=b)))
            }

        if (parms$case==2) {
            names(newtheta) <- 'Intercept'
            names(b) <- parms$levellist[[1]]
            theta <- list(newtheta)
            names(theta) <- parms$gname
            b <- list(b)
            names(b) <- parms$gname
            return(list(theta=theta, b=b))
            }
        if (parms$case==3) {
            ngroup <- length(parms$levellist)
            theta <- vector('list', ngroup)
            names(theta) <- parms$gname
            for (i in 1:ngroup) 
                theta[[parms$gname[i]]] <- c('(Intercept)'=newtheta[i])

            if (parms$collapse) {
                names(b) <- parms$levellist[[1]]
                random <- list(b)
                names(random) <- parms$gname[[1]]
                }
            else {
                names(b) <- unlist(parms$levellist)
                random <- split(b, rep(1:ngroup, parms$nlevel))
                names(random) <- parms$gname
                }
            return(list(theta=theta, b=random))
            }
        if (parms$case==4) {
            intercept <- parms$intercept
            ngroup <- length(parms$nlevel)
            nvar <- parms$nvar

            # Deal with b
            random <- split(b, rep(rep(1:ngroup, parms$nlevel), intercept +nvar))
            names(random) <- parms$gname
            if (intercept) {
                colname <- c("Intercept", parms$xname)
                }
            else {
                colname <- parms$xname
                }

            for (i in 1:ngroup) {
                temp<- matrix(random[[i]], ncol=length(colname))
                random[[i]] <- temp 
                dimnames(random[[i]]) <- list(parms$levellist[[i]], colname)
                }
            
            # Deal with theta
            tfun <- function(x, n= 1 + nvar) {
                tmat <- matrix(0., n, n)
                tmat[row(tmat) >= col(tmat)] <- x
                offdiag <- row(tmat) > col(tmat)
                tmat[offdiag] <- (tmat[offdiag]-1)/(tmat[offdiag]+1)
                dimnames(tmat) <- list(colname, colname)
                tmat + t(tmat) - diag(diag(tmat))
                }
            parms.per.group <- length(newtheta)/ngroup
            if (parms.per.group==1) { #nvar=1, intercept=F case
                theta <- as.list(newtheta)
                theta <- lapply(theta, function(x) {names(x)<- parms$xname; x})
                names(theta) <- parms$gname
                }
            else {
                theta <- vector('list', ngroup)
                names(theta) <- parms$gname
                for (i in 1:ngroup) 
                    theta[[i]] <- tfun(newtheta[1:parms.per.group + 
                                             parms.per.group*(i-1)])
                }
            return(list(theta=theta, b=random))
            }
        }
    out <- list(initialize=initialize, generate=generate, wrapup=wrapup)
    oldClass(out) <- 'coxmevar'
    out
    }
