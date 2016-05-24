# Automatically generated from all.nw using noweb
coxme <- function(formula,  data, 
        weights, subset, na.action, init, 
        control, ties= c("efron", "breslow"),
        varlist, vfixed, vinit,
        x=FALSE, y=TRUE, 
        refine.n=0, random, fixed, variance,  ...) {

    #time0 <- proc.time()    #debugging line
    ties <- match.arg(ties)
    Call <- match.call()

    if (!missing(fixed)) {
        if (missing(formula)) {
            formula <- fixed
            warning("The 'fixed' argument of coxme is depreciated")
            }
        else stop("Both a fixed and a formula argument are present")
        }
    if (!missing(random)) {
        warning("The random argument of coxme is depreciated")
        if (class(random) != 'formula' || length(random) !=2) 
            stop("Invalid random formula")
        j <- length(formula)   #will be 2 or 3, depending on if there is a y

        # Add parens to the random formula and paste it on
        formula[[j]] <- call('+', formula[[j]], call('(', random[[2]]))  
        }

    if (!missing(variance)) {
        warning("The variance argument of coxme is depreciated")
        vfixed <- variance
        }

    temp <- call('model.frame', formula= subbar(formula))
    for (i in c('data', 'subset', 'weights', 'na.action'))
        if (!is.null(Call[[i]])) temp[[i]] <- Call[[i]]
    if (is.R()) m <- eval.parent(temp)
    else        m <- eval(temp, sys.parent())
        Y <- model.extract(m, "response")
        n <- nrow(Y)
        if (n==0) stop("No observations remain in the data set")
        if (!inherits(Y, "Surv")) stop("Response must be a survival object")
        type <- attr(Y, "type")
        if (type!='right' && type!='counting')
            stop(paste("Cox model doesn't support '", type,
                              "' survival data", sep=''))

        weights <- model.weights(m)
        if (length(weights) ==0) weights <- rep(1.0, n)
        else if (any(weights <=0))
            stop("Negative or zero weights are not allowed")

        offset <- model.offset(m)
        if (length(offset)==0) offset <- rep(0., n)

        # Check for penalized terms; the most likely is pspline
        pterms <- sapply(m, inherits, 'coxph.penalty')
        if (any(pterms)) {
            stop("You cannot have penalized terms in coxme")
            }

        if (missing(control)) control <- coxme.control(...)
        if (missing(init)) init <- NULL
        flist <- formula1(formula)
        if (hasAbar(flist$fixed))
            stop("Invalid formula: a '|' outside of a valid random effects term")

        special <- c("strata", "cluster")
        Terms <- terms(flist$fixed, special)
        attr(Terms,"intercept")<- 1  #Cox model always has \Lambda_0
        strats <- attr(Terms, "specials")$strata
        cluster<- attr(Terms, "specials")$cluster
        if (length(cluster)) {
            stop ("A cluster() statement is invalid in coxme")
            }
        if (length(strats)) {
            temp <- untangle.specials(Terms, 'strata', 1)
            if (length(temp$vars)==1) strata.keep <- m[[temp$vars]]
            else strata.keep <- strata(m[,temp$vars], shortlabel=T)
            strats <- as.numeric(strata.keep)
            X <- model.matrix(Terms[-temp$terms], m)[,-1,drop=F]
            }
        else X <- model.matrix(Terms, m)[,-1,drop=F]


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
    fit <- coxme.fit(X, Y, strats, offset, init, control, weights=weights,
                     ties=ties, row.names(m),
                     fmat, zmat, varlist, vparm, 
                     itheta, ntheta, ncoef, refine.n,
                     is.variance = any(is.variance))
    if (is.character(fit)) {
        fit <- list(fail=fit)
        oldClass(fit) <- 'coxme'
        return(fit)
        }
    fcoef <- fit$coefficients
    nvar <- length(fcoef)
    if (length(fcoef)>0 && any(is.na(fcoef))) {
        vars <- (1:length(fcoef))[is.na(fcoef)]
        msg <-paste("X matrix deemed to be singular; variable",
                        paste(vars, collapse=" "))
        warning(msg)
        }
    if (length(fcoef) >0) {
        names(fit$coefficients) <- dimnames(X)[[2]]
        }
    rlinear <- rep(0., nrow(Y))
    if (length(fmat)) {
        for (i in 1:ncol(fmat)) {
            rlinear <- rlinear + fit$frail[fmat[,i]]
            }
        }
    if (length(zmat)) {
        indx <- if (length(fmat)>0) max(fmat) else 0
        for (i in 1:ncol(zmat))
            rlinear <- rlinear + fit$frail[indx+i]*zmat[,i]
        }

    if (nvar==0) fit$linear.predictor <- rlinear
    else fit$linear.predictor <- as.vector(rlinear + c(X %*% fit$coef))
    newtheta <- random.coef <- list()  
    nrandom <- length(varlist)
    sindex <- rep(1:nrandom, ntheta) # which thetas to which terms
    bindex <- rep(row(ncoef), ncoef) # which b's to which terms
    for (i in 1:nrandom) {
        temp <- varlist[[i]]$wrapup(fit$theta[sindex==i], fit$frail[bindex==i], 
                                    vparm[[i]])
        newtheta <- c(newtheta, temp$theta)
        if (!is.list(temp$b)) {
            temp$b <- list(temp$b)
            names(temp$b) <- paste("Random", i, sep='')
            }
        random.coef <- c(random.coef, temp$b)
        }
    fit$frail <- random.coef

    fit$vcoef <- newtheta
    fit$theta <- NULL
    fit$n <- c(sum(Y[,ncol(Y)]), nrow(Y))
    fit$terms <- Terms
    fit$assign <- attr(X, 'assign')
    fit$formulaList <- flist

    na.action <- attr(m, "na.action")
    if (length(na.action)) fit$na.action <- na.action
    if (x)  {
        fit$x <- X
        if (length(strats)) fit$strata <- strata.keep
        }
    if (y)     fit$y <- Y
    if (!is.null(weights) && any(weights!=1)) fit$weights <- weights

    fit$formula <- as.vector(attr(Terms, "formula"))
    fit$call <- Call
    fit$ties <- ties
    names(fit$loglik) <- c("NULL", "Integrated", "Penalized")
    oldClass(fit) <- 'coxme'
    fit
}
