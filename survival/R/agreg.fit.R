# Automatically generated from the noweb directory
agreg.fit <- function(x, y, strata, offset, init, control,
                        weights, method, rownames)
    {
    n <- nrow(y)
    nvar <- ncol(x)
    start <- y[,1]
    stopp <- y[,2]
    event <- y[,3]
    if (all(event==0)) stop("Can't fit a Cox model with 0 failures")

    # Sort the data (or rather, get a list of sorted indices)
    #  For both stop and start times, the indices go from last to first
    if (length(strata)==0) {
        sort.end  <- order(-stopp) -1L #indices start at 0 for C code
        sort.start<- order(-start) -1L
        newstrat  <- n
        }
    else {
        sort.end  <- order(strata, -stopp) -1L
        sort.start<- order(strata, -start) -1L
        newstrat  <- cumsum(table(strata))
        }
    if (missing(offset) || is.null(offset)) offset <- rep(0.0, n)
    if (missing(weights)|| is.null(weights))weights<- rep(1.0, n)
    else if (any(weights<=0)) stop("Invalid weights, must be >0")
    else weights <- as.vector(weights)

    if (is.null(nvar) || nvar==0) {
        # A special case: Null model.  Just return obvious stuff
        #  To keep the C code to a small set, we call the usual routines, but
        #  with a dummy X matrix and 0 iterations
        nvar <- 1
        x <- matrix(as.double(1:n), ncol=1)  #keep the .C call happy
        maxiter <- 0
        nullmodel <- TRUE
        if (length(init) !=0) stop("Wrong length for inital values")
        init <- 0.0  #dummy value to keep a .C call happy (doesn't like 0 length)
        }
    else {
        nullmodel <- FALSE
        maxiter <- control$iter.max
        
        if (is.null(init)) init <- rep(0., nvar)
        if (length(init) != nvar) stop("Wrong length for inital values")
        }

    # the returned value of agfit$coef starts as a copy of init, so make sure
    #  is is a vector and not a matrix; as.double does so.
    # Solidify the storage mode of other arguments
    storage.mode(y) <- storage.mode(x) <- "double"
    storage.mode(offset) <- storage.mode(weights) <- "double"
    storage.mode(newstrat) <- "integer"
    agfit <- .Call(Cagfit4, 
                   y, x, newstrat, weights, 
                   offset,
                   as.double(init), 
                   sort.start, sort.end, 
                   as.integer(method=="efron"),
                   as.integer(maxiter), 
                   as.double(control$eps),
                   as.double(control$toler.chol),
                   as.integer(1)) # internally rescale

    var <- matrix(agfit$imat,nvar,nvar)
    coef <- agfit$coef
    if (agfit$flag[1] < nvar) which.sing <- diag(var)==0
    else which.sing <- rep(FALSE,nvar)

    infs <- abs(agfit$u %*% var)
    if (maxiter >1) {
        if (agfit$iter > maxiter)
            warning("Ran out of iterations and did not converge")
        else {
            infs <- ((infs > control$eps) & 
                     infs > control$toler.inf*abs(coef))
            if (any(infs))
                warning(paste("Loglik converged before variable ",
                              paste((1:nvar)[infs],collapse=","),
                                          "; beta may be infinite. "))
        }
    }
    lp  <- as.vector(x %*% coef + offset - sum(coef * colMeans(x)))
    score <- as.double(exp(lp))
    resid <- .Call(Cagmart3,
                   y, score, weights,
                   newstrat,
                   cbind(sort.end, sort.start),
                   as.integer(method=='efron'))
    names(resid) <- rownames
    if (nullmodel) {
        list(loglik=agfit$loglik[2],
             linear.predictors = offset,
             residuals = resid,
             method= c("coxph.null", 'coxph') )
    }
    else {
        names(coef) <- dimnames(x)[[2]]
        if (maxiter > 0) coef[which.sing] <- NA  # always leave iter=0 alone
        flag <- agfit$flag
        names(flag) <- c("rank", "rescale", "step halving")
        
        concordance <- survConcordance.fit(y, lp, strata, weights) 
        list(coefficients  = coef,
             var    = var,
             loglik = agfit$loglik,
             score  = agfit$sctest,
             iter   = agfit$iter,
             linear.predictors = as.vector(lp),
             residuals = resid,
             means = colMeans(x),
             concordance = concordance,
             first = agfit$u,
             info = flag,
             method= 'coxph')
    }
}  
