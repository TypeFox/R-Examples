# Automatically generated from all.nw using noweb
coxmeMlist <- function(varlist, rescale=FALSE, pdcheck=TRUE,  positive=TRUE) {
    # Because of environments, the init function will inherit the
    #  four variables below 
    varlist <- varlist
    rescale <- rescale
    pdcheck <- pdcheck
    if (!is.logical(positive)) stop("Invalid value for postive argument")
    positive <- positive
    initialize <- function(vinit, fixed, intercept, G, X, control) {
        ngroup <- min(length(G), ncol(G))
        nvar   <- min(length(X), ncol(X))  # a NULL or a nx0 matrix yields 0
        if (ngroup >0 & nvar >0)
            return(list(error="Mlist cannot have both covariates and grouping"))
        if (!is.list(varlist)) varlist <- list(varlist)  # a naked matrix
        noname <- all(sapply(varlist, function(x) is.null(dimnames(x)) || 
                          (is.null(dimnames(x)[[1]]) & is.null(dimnames(x)[[2]]))))
        namefun <- function(x, names) {
            if (all(dim(x)== rep(length(names),2))) 
                dimnames(x) <- list(names, names)
            x
            }
        if (ngroup >0) {
            n <- nrow(G)
            G <- expand.nested(G)
            groups <- G[[ngroup]]  #drop all but the last
            if (noname) varlist <- lapply(varlist, namefun, levels(groups))
            if (any(sapply(varlist, function(x) inherits(x, "Matrix"))))
                varlist <- lapply(varlist, function(x) as(x, "bdsmatrix"))
            tlist <- bdsmatrix.reconcile(varlist, levels(groups))
            bname <-  dimnames(tlist[[1]])[[1]]
            imap <- matrix(match(groups, bname))
            xmap <- NULL
            rname <- names(G)[[ngroup]]
            }
        else {
            n <- nrow(X)
            bname <- dimnames(X)[[2]]
            if (noname) varlist <- lapply(varlist, namefun, bname)
            tlist <- bdsmatrix.reconcile(varlist, bname)
            # sparse matrices (bdsmatrix or Matrix) are illegal, for now, 
            #   for covariates
            tlist <- lapply(tlist, as.matrix)
            xmap <- match(dimnames(X)[[2]], bname)
            xmap <- matrix(rep(xmap, n), nrow=n, byrow=T)
            imap <- NULL
            rname <- "(Shrink)"
            }

        ntheta <- length(varlist)
        fudge <- seq(1, 1.5, length=ntheta)
        itheta <- vector('list', ntheta)
        for (i in 1:ntheta) itheta[[i]] <- control$varinit * fudge[i]

        if (length(vinit) >0) {
            if (length(vinit) != ntheta)
                return(list(error="Wrong length for initial values"))
            indx <- !is.na(vinit) & vinit !=0  #which to use
            if (any(indx)) itheta[indx] <- vinit[indx]
            }

        which.fixed <- rep(FALSE, ntheta)
        if (length(fixed) >0) {
            if (length(fixed) != ntheta)
                return(list(error="Wrong length for fixed values"))
            indx <- !is.na(fixed) & fixed !=0  #which to use
            if (any(indx)) {
                itheta[indx] <- fixed[indx]
                which.fixed[indx] <- TRUE
                }
            }

        if (length(positive)==1) positive <- rep(positive, ntheta)
        if (length(positive) != ntheta)
            return(list(error="Wrong length for positive parameter"))
        if (any(unlist(itheta[positive]) <=0))
            return(list(error="Invalid initial value, must be positive"))        
        itheta[positive] <- lapply(itheta[positive], log)    
        for (j in 1:ntheta) {
            kmat <- tlist[[j]]
            if (rescale) {
                temp <- diag(kmat)
                if (any(temp==0))
                    return(list(error="Diagonal of a variance matrix is zero"))
        #        if (any(temp != temp[1])) 
        #            warning("Diagonal of variance matrix is not constant")
                if (max(temp) !=1) {
                    kmat <- kmat/max(temp)
                    tlist[[j]] <- kmat
                    }
                }
            if (pdcheck) {
                temp <- gchol(kmat)
                if (any(diag(temp) < 0))
                    return(list(error="A variance matrix is not non-negative definite"))
                }
            }

        # itheta is a list with vectors of initial values
        # theta is a vector, and only the fixed values need to be correct (the others
        #  are replaced by the parent routine).  All fixed "inits" are of length 1.
        theta <- sapply(itheta, function(x) x[1])    
        list(theta=itheta[!which.fixed], imap=imap, X=X, xmap=xmap,
             parms=list(varlist=tlist, theta=theta, fixed=which.fixed,
                        bname=bname, rname=rname, positive=positive,
                        vname=names(varlist)))
        }
     generate <- function(newtheta, parms) {
         theta <- parms$theta
         theta[!parms$fixed] <- newtheta
         if (any(parms$positive)) theta[parms$positive] <- 
              exp(pmax(-36, pmin(36, theta[parms$positive])))

         varmat <- parms$varlist[[1]] * theta[1]
         if (length(theta) >1) {
             for (i in 2:length(theta)) {
                 varmat <- varmat + theta[i]*parms$varlist[[i]]
                 }
             }
         varmat
         }
    wrapup <- function(newtheta, b, parms) {
            theta <- parms$theta
            theta[!parms$fixed] <- newtheta
            theta[parms$positive] <- exp(theta[parms$positive])

            defaultname <- paste("Vmat", 1:length(theta), sep=".")
            vname <- parms$vname
            if (length(vname)==0) vname <- defaultname
            else if (any(vname=='')){
                indx <- which(vname=='')
                vname[indx] <- defaultname[indx]
                }
            names(theta) <- vname
            theta <- list(theta)
            names(theta) <- parms$rname
            names(b) <- parms$bname
            b <- list(b)
            names(b) <- parms$rname
            list(theta=theta, b=b)
            }
    out <- list(initialize=initialize, generate=generate, wrapup=wrapup)
    class(out) <- 'coxmevar'
    out
    }
