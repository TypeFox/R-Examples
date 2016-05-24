# Automatically generated from all.nw using noweb
lmekin <- function(formula,  data, 
        weights, subset, na.action,
        control, varlist, vfixed, vinit, 
        method=c("ML", "REML"),
        x=FALSE, y=FALSE, model=FALSE,
        random, fixed, variance,  ...) {

    Call <- match.call()
    sparse <- c(1,0)  #needed for compatablily with coxme code
    if (!missing(fixed)) {
        if (missing(formula)) {
            formula <- fixed
            warning("The 'fixed' argument of lmekin is depreciated")
            }
        else stop("Both a fixed and a formula argument are present")
        }
    if (!missing(random)) {
        warning("The random argument of lmekin is depreciated")
        if (class(random) != 'formula' || length(random) !=2) 
            stop("Invalid random formula")
        j <- length(formula)   #will be 2 or 3, depending on if there is a y

        # Add parens to the random formula and paste it on
        formula[[j]] <- call('+', formula[[j]], call('(', random[[2]]))  
        }

    if (!missing(variance)) {
        warning("The variance argument of lmekin is depreciated")
        vfixed <- variance
        }

    method <- match.arg(method)

    temp <- call('model.frame', formula= subbar(formula))
    for (i in c('data', 'subset', 'weights', 'na.action'))
        if (!is.null(Call[[i]])) temp[[i]] <- Call[[i]]
    m <- eval.parent(temp)

    Y <- model.extract(m, "response")
    n <- length(Y)
    if (n==0) stop("data has no observations")

    weights <- model.weights(m)
    if (length(weights) ==0) weights <- rep(1.0, n)
    else if (any(weights <=0))
        stop("Negative or zero weights are not allowed")

    offset <- model.offset(m)
    if (length(offset)==0) offset <- rep(0., n)

    # Check for penalized terms; the most likely is pspline
    pterms <- sapply(m, inherits, 'coxph.penalty')
    if (any(pterms)) {
        stop("You cannot have penalized terms in lmekin")
        }

    if (missing(control)) control <- lmekin.control(...)
    flist <- formula1(formula)
    if (hasAbar(flist$fixed))
        stop("Invalid formula: a '|' outside of a valid random effects term")

    special <- c("strata", "cluster")
    Terms <- terms(flist$fixed, special)
    if (length(attr(Terms, "specials")$strata))
        stop ("A strata term is invalid in lmekin")
    if (length(attr(Terms, "specials")$cluster))
        stop ("A cluster term is invalid in lmekin")
    X <- model.matrix(Terms, m)
    nrandom <- length(flist$random)
    if (nrandom ==0) stop("No random effects terms found")
    vparm <- vector('list', nrandom)
    is.variance <- rep(TRUE, nrandom)  #penalty fcn returns a variance or penalty?
    ismat <- function (x) {
        inherits(x, "matrix") || inherits(x, "bdsmatrix") | inherits(x, "Matrix")
    }
    if (missing(varlist) || is.null(varlist)) {
        varlist <- vector('list', nrandom)
        for (i in 1:nrandom) varlist[[i]] <- coxmeFull() #default
        }
    else {
        if (is.function(varlist)) varlist <- varlist()
        if (class(varlist)=='coxmevar') varlist <- list(varlist)
        else if (ismat(varlist))
            varlist <- list(coxmeMlist(list(varlist)))
        else {
            if (!is.list(varlist)) stop("Invalid varlist argument")
            if (all(sapply(varlist, ismat))) {
                # A list of matrices
                if (nrandom >1) 
                    stop(paste("An unlabeled list of matrices is",
                               "ambiguous when there are multiple random terms"))
                else varlist <- list(coxmeMlist(varlist))
                }
            else {  #the user gave me a list, not all matrices
                for (i in 1:length(varlist)) {
                    if (is.function(varlist[[i]])) 
                        varlist[[i]] <-varlist[[i]]()
                    if (ismat(varlist[[i]]))
                        varlist[[i]] <- coxmeMlist(list(varlist[[i]]))
                    if (class(varlist[[i]]) != 'coxmevar') {
                        if (is.list(varlist[[i]])) {
                            if (all(sapply(varlist[[i]], ismat)))
                                varlist[[i]] <- coxmeMlist(varlist[[i]])
                            else stop("Invalid varlist element")
                            }
                        else stop("Invalid varlist element")
                        }
                    }
                }
            }
        while(length(varlist) < nrandom) varlist <- c(varlist, list(coxmeFull()))
        }


    if (!is.null(names(varlist))) { # put it in the right order
        vname <- names(varlist)
        stop("Cannot (yet) have a names varlist")
        indx <- pmatch(vname, names(random), nomatch=0)
        if (any(indx==0 & vname!=''))
            stop(paste("Varlist element not matched:", vname[indx==0 & vname!='']))
        if (any(indx>0)) {
            temp <- vector('list', nrandom)
            temp[indx] <- varlist[indx>0]
            temp[-indx]<- varlist[indx==0]
            varlist <- temp
            }
        }
        
    #check validity (result is never used)
    check <- sapply(varlist, function(x) {
           fname <- c("initialize", "generate", "wrapup")
           indx <- match(fname, names(x))
           if (any(is.na(x)))
               stop(paste("Member not found in variance function:",
                          fname(is.na(indx))))
           if (length(x) !=3 || any(!sapply(x, is.function)))
               stop("Varlist objects must consist of exaclty three functions")
       })

    getcmat <- function(x, mf) {
        if (is.null(x) || x==1) return(NULL)
        Terms <- terms(eval(call("~", x)))
        attr(Terms, 'intercept') <- 0  #ignore any "1+" that is present

        varnames <-  attr(Terms, 'term.labels')
        ftemp <- sapply(mf[varnames], is.factor)
        if (any(ftemp)) {
            clist <- lapply(mf[varnames[ftemp]], 
                            function(x) diag(length(levels(x))))
            model.matrix(Terms, mf, contrasts.arg =clist)
            }
        else model.matrix(Terms, mf)
        }
    getGroupNames <- function(x) {
        if (is.call(x) && x[[1]]==as.name('/')) 
            c(getGroupNames(x[[2]]), getGroupNames(x[[3]]))
        else deparse(x)
        }

    getgroups <- function(x, mf) {
        if (is.numeric(x)) return(NULL)  # a shrinkage effect like (x1+x2 | 1)
        varname <- getGroupNames(x)
        indx <- match(varname, names(mf), nomatch=0)
        if (any(indx==0)) stop(paste("Invalid grouping factor", varname[indx==0]))
        else data.frame(lapply(mf[indx], as.factor))
        }
    if (missing(vinit) || is.null(vinit)) vinit <- vector('list', nrandom)
    else {
        if (nrandom==1) {
            if (is.numeric(vinit)) vinit <- list(vinit)
            else if (is.list(vinit)) vinit <- list(unlist(vinit))
        }
        if (!is.list(vinit)) stop("Invalid value for `vinit` parameter")
        if (length(vinit) > nrandom) 
            stop (paste("Vinit must be a list of length", nrandom))
        if (!all(sapply(vinit, function(x) (is.null(x) || is.numeric(x))))) 
            stop("Vinit must contain numeric values") 
        
        if (length(vinit) < nrandom) 
            vinit <- c(vinit, vector('list', nrandom - length(vinit)))
                       
        tname <- names(vinit)
        if (!is.null(tname)) {
            stop("Named initial values not yet supported")
            #temp <- pmatch(tname, names(flist$random), nomatch=0)
            #temp <- c(temp, (1:nrandom)[-temp])
            #vinit <- vinit[temp]
            }
      }

    if (missing(vfixed) || is.null(vfixed)) vfixed <- vector('list', nrandom)
    else {
        if (nrandom==1) {
            if (is.numeric(vfixed)) vfixed <- list(vfixed)
            else if (is.list(vfixed)) vfixed <- list(unlist(vfixed))
        }
        if (!is.list(vfixed)) stop("Invalid value for `vfixed` parameter")
        if (length(vfixed) > nrandom) 
            stop (paste("Vfixed must be a list of length", nrandom))
        if (!all(sapply(vfixed,  function(x) (is.null(x) || is.numeric(x))))) 
            stop("Vfixed must contain numeric values") 

        if (length(vfixed) < nrandom) 
            vfixed <- c(vfixed, vector('list', nrandom - length(vfixed)))
                       
        tname <- names(vfixed)
        if (!is.null(tname)) {
            temp <- pmatch(tname, names(flist$random), nomatch=0)
            temp <- c(temp, (1:nrandom)[-temp])
            vfixed <- vfixed[temp]
            }
      }
    newzmat <- function(xmat, xmap) {
        n <- nrow(xmap)
        newz <- matrix(0., nrow=n, ncol=max(xmap))
        for (i in 1:ncol(xmap)) 
            newz[cbind(1:n, xmap[,i])] <- xmat[,i]
        newz
        }
    fmat <- zmat <- matrix(0, nrow=n, ncol=0)
    ntheta <- integer(nrandom)
    ncoef  <- matrix(0L, nrandom, 2, dimnames=list(NULL, c("intercept", "slope")))
    itheta <-  NULL   #initial values of parameters to iterate over

    for (i in 1:nrandom) {
        f2 <- formula2(flist$random[[i]])
        if (f2$intercept & f2$group==1)
            stop(paste("Error in random term ", i, 
                       ": Random intercepts require a grouping variable", sep=''))
        vfun <- varlist[[i]]
        if (!is.null(f2$interaction)) stop("Interactions not yet written")

        cmat <- getcmat(f2$fixed, m)
        groups <- getgroups(f2$group, m)
        ifun <- vfun$initialize(vinit[[i]], vfixed[[i]], intercept=f2$intercept, 
                            groups, cmat, control)
        if (!is.null(ifun$error)) 
            stop(paste("In random term ", i, ": ", ifun$error, sep=''))
        vparm[[i]] <- ifun$parms
        if (!is.null(ifun$is.variance)) is.variance[i] <- ifun$is.variance
        itheta <- c(itheta, ifun$theta)
        ntheta[i] <- length(ifun$theta)

        if (f2$intercept) {
            if (!is.matrix(ifun$imap) || nrow(ifun$imap) !=n) 
                stop(paste("In random term ", i, 
                           ": Invalid intercept matrix F", sep=''))
            temp <- sort(unique(c(ifun$imap)))
            if (any(temp != 1:length(temp)))
                stop(paste("In random term ", i,
                           ": intercept matrix has an invalid element", sep=''))

            if (ncol(fmat) >0) fmat <- cbind(fmat, ifun$imap + max(fmat))
            else fmat <- ifun$imap
            ncoef[i,1] <- 1+ max(ifun$imap) - min(ifun$imap)
            }

        if (length(cmat)>0) {
            if (is.null(ifun$xmap) || is.null(ifun$X) ||
                !is.matrix(ifun$xmap) || !is.matrix(ifun$X) ||
                nrow(ifun$xmap) !=n || nrow(ifun$X) != n ||
                ncol(ifun$xmap) != ncol(ifun$X))
                stop(paste("In random term ", i,
                           "invalid X/xmap pair"))
            if (f2$intercept) xmap <- ifun$xmap - max(ifun$imap)
            else xmap <- ifun$xmap
            if (any(sort(unique(c(xmap))) != 1:max(xmap)))
                 stop(paste("In random term ", i,
                           ": xmap matrix has an invalid element", sep=''))
            
            temp <- newzmat(ifun$X, xmap)
            ncoef[i,2] <- ncol(temp)
            zmat <- cbind(zmat, temp)
            }
    }
    if (any(is.variance) & !all(is.variance))
             stop("All variance terms must have the same is.variance setting") 
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
    #Define Z^* and X^*
    itemp <- split(row(fmat), fmat)
    zstar1 <- new("dgCMatrix", 
                  i= as.integer(unlist(itemp) -1),
                  p= as.integer(c(0, cumsum(unlist(lapply(itemp, length))))),
                  Dim=as.integer(c(n, max(fmat))),
                  Dimnames= list(NULL, NULL),
                  x= rep(1.0, length(fmat)),
                  factors=list())
    if (length(zmat) >0)  {
        # there were random slopes as well
        zstar1 <- cBind(zstar1, as(Matrix(zmat), "dgCMatrix"))
    }

    nfrail <- ncol(zstar1)
    nvar <- ncol(X)
    if (nvar == 0)  xstar <- NULL  #model with no covariates
    else xstar <- rBind(Matrix(X, sparse=FALSE),
                        matrix(0., nrow=nfrail, ncol=ncol(X)))
    ystar <- c(Y, rep(0.0, nfrail))
    mydiag <- function(x) {
        if (class(x)=="sparseQR") diag(x@R)
        else diag(qr.R(x))
    }

    logfun <- function(theta, best=0) {
        vmat <- kfun(theta, varlist, vparm, ntheta, ncoef)
        Delta <- t(solve(chol(as(vmat, "dsCMatrix"), pivot=FALSE)))
        zstar <- rBind(zstar1, Delta)
        qr1 <- qr(zstar)
        dd <- mydiag(qr1) 
        cvec <- as.vector(qr.qty(qr1, ystar))[-(1:nfrail)]  #residual part
        if (nvar >0) {  # have covariates
            qr2 <- qr(qr.qty(qr1, xstar)[-(1:nfrail),])
            cvec <- qr.qty(qr2, cvec)[-(1:nvar)]  #residual part
            if (method!= "ML") dd <- c(dd, mydiag(qr2))
        }

        loglik <- sum(log(abs(diag(Delta)))) - sum(log(abs(dd)))
        if (method=="ML") loglik <- loglik - .5*n*log(sum(cvec^2))
        else              loglik <- loglik - .5*length(cvec)*log(sum(cvec^2))
        
        best - (loglik+1)  #optim() wants to minimize rather than maximize
    }

    nstart <- sapply(itheta, length)
    if (length(nstart) ==0) theta <- NULL #no thetas to solve for
    else {
        #iteration is required
        #make a matrix of all possible starting estimtes
        testvals <- do.call(expand.grid, itheta)
        bestlog <- NULL
        for (i in 1:nrow(testvals)) {
            ll <- logfun(as.numeric(testvals[i,]))
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

        optpar <- control$optpar
        optpar$hessian <- TRUE
        mfit <- do.call('optim', c(list(par= theta, fn=logfun, gr=NULL,
                                        best=bestlog), optpar))
        theta <- mfit$par
    }
    vmat <-  kfun(theta, varlist, vparm, ntheta, ncoef)
    Delta <- t(solve(chol(as(vmat, "dsCMatrix"), pivot=FALSE)))
    zstar <- rBind(zstar1, Delta)
    qr1 <- qr(zstar)
    dd <- mydiag(qr1)
    ctemp <- as.vector(qr.qty(qr1, ystar))
    cvec <- ctemp[-(1:nfrail)]  #residual part

    if (is.null(xstar)) { #No X covariates
        rcoef <- qr.coef(qr1, ystar)
        yhat <- qr.fitted(qr1, ystar)
    }
    else {
        qtx <- qr.qty(qr1, xstar)
        qr2 <- qr(qtx[-(1:nfrail),,drop=F])
        if (method!="ML") dd <- c(dd, mydiag(qr2))

        fcoef <-qr.coef(qr2, cvec)
        yresid <- ystar - xstar %*% fcoef
        rcoef <- qr.coef(qr1, yresid)
        cvec <- qr.qty(qr2, cvec)[-(1:nvar)] #residual part
        if (class(qr2)=="sparseQR") varmat <- chol2inv(qr2@R)
        else varmat <- chol2inv(qr.R(qr2))
        yhat <- as.vector(zstar1 %*% rcoef + X %*% fcoef) #kill any names
    }

    if (method=="ML") {
        sigma2 <- sum(cvec^2)/n  #MLE estimate  
        loglik <- sum(log(abs(diag(Delta)))) - 
              (sum(log(abs(dd))) + .5*n*(log(2*pi) +1 + log(sigma2)))
    }
    else {
        np <- length(cvec)  # n-p
        sigma2 <- mean(cvec^2)  # divide by n-p
        loglik <- sum(log(abs(diag(Delta)))) -
            (sum(log(abs(dd))) + .5*np*(log(2*pi) + 1+ log(sigma2)))
    }
            
    # Debugging code, set the argument to TRUE only during testing
    if (FALSE) {
        # Compute the alternate way (assumes limited reordering)
        zx <- cBind(zstar, as(xstar, class(zstar)))
        qr3 <- qr(zx)
        cvec3 <- qr.qty(qr3, ystar)[-(1:(nvar+nfrail))]
        if (method=="ML")  dd3 <- (diag(myqrr(qr3)))[1:nfrail]
        else               dd3 <- (diag(myqrr(qr3)))[1:(nfrail+nvar)]
        #all.equal(dd, dd3)
        #all.equal(cvec, cvec3)
        acoef <- qr.coef(qr3, ystar)
        browser()
    }

    newtheta <- random.coef <- list()  
    nrandom <- length(varlist)
    sindex <- rep(1:nrandom, ntheta) #which thetas to which terms
    bindex <- rep(1:nrandom, rowSums(ncoef)) # which b's to which terms
    for (i in 1:nrandom) {
        temp <- varlist[[i]]$wrapup(theta[sindex==i], rcoef[bindex==i], 
                                    vparm[[i]])
        newtheta <- c(newtheta, lapply(temp$theta, function(x) x*sigma2))
        if (!is.list(temp$b)) {
            temp$b <- list(temp$b)
            names(temp$b) <- paste("Random", i, sep='')
            }
        random.coef <- c(random.coef, temp$b)
        }

    if (length(fcoef) >0) {
        # There are fixed effects
        nvar <- length(fcoef)
        fit <- list (coefficients=list(fixed=fcoef, random=random.coef),
                     var = varmat * sigma2,
                     vcoef =newtheta,
                     residuals= Y- yhat,
                     method=method,
                     loglik=loglik,
                     sigma=sqrt(sigma2),
                     n=n,
                     call=Call)
    }
    else fit <- list(coefficients=list(fixed=NULL, random=random.coef),
                     vcoef=newtheta,
                     residuals=Y - yhat,
                     method=method,
                     loglik=loglik,
                     sigma=sqrt(sigma2),
                     n=n,
                     call=Call)

    if (!is.null(theta)) {
        fit$rvar <- mfit$hessian
        fit$iter <- mfit$counts
    }
    if (x) fit$x <- X
    if (y) fit$y <- Y
    if (model) fit$model <- m

    na.action <- attr(m, "na.action")
    if (length(na.action)) fit$na.action <- na.action
                             
    class(fit) <- "lmekin"
    fit
    }
residuals.lmekin <- function(object, ...) {
    if (length(object$na.action)) naresid(object$.na.action, object$residuals)
    else object$residuals
}
