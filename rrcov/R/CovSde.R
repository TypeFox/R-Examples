CovSde <- function(x,
                   nsamp,
                   maxres,
                   tune=0.95,
                   eps=0.5,
                   prob=0.99,
                   seed=NULL,
                   trace=FALSE,
                   control)
{

    ## Analize and validate the input parameters ...

    ## if a control object was supplied, take the option parameters from it,
    ## but if single parameters were passed (not defaults) they will override the
    ## control object.

    if(!missing(control)){
        defcontrol <- CovControlSde()       # default control
        ##  no default if(nsamp == defcontrol@nsamp)       nsamp <- control@nsamp
        ##  no default if(maxres == defcontrol@maxres)     maxres <- control@maxres
        if(tune == defcontrol@tune)         tune <- control@tune
        if(eps == defcontrol@eps)           eps <- control@eps
        if(prob == defcontrol@prob)         prob <- control@prob
        if(is.null(seed) || seed == defcontrol@seed)         seed <- control@seed
        if(trace == defcontrol@trace)       trace <- control@trace
    }

    if(length(seed) > 0) {
        if(exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE))  {
            seed.keep <- get(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
            on.exit(assign(".Random.seed", seed.keep, envir=.GlobalEnv))
        }
        assign(".Random.seed", seed, envir=.GlobalEnv)
    }

    if(is.data.frame(x))
        x <- data.matrix(x)
    else if (!is.matrix(x))
        x <- matrix(x, length(x), 1,
         dimnames = list(names(x), deparse(substitute(x))))
    xcall <- match.call()

    ## drop all rows with missing values (!!) :
    na.x <- !is.finite(x %*% rep(1, ncol(x)))
    ok <- !na.x
    x <- x[ok, , drop = FALSE]
    dx <- dim(x)
    if(!length(dx))
        stop("All observations have missing values!")
    dimn <- dimnames(x)
    n <- dx[1]
    p <- dx[2]

    if(missing(nsamp))
        nsamp <- ceiling(log(1 - prob)/log(1 - (1 - eps)^(p+1)))
    else if(!is.numeric(nsamp))
        stop("If present, nsamp must be a nonnegative integer or")

    if(nsamp != 0)
        nsamp <- max(1000, nsamp)
    nsamp <- min(nsamp, .Machine$integer.max)

    if(missing(maxres))
        maxres <- 2 * nsamp
    else if(!is.numeric(maxres))
        stop(sQuote("maxres"), " is not a positive integer")
    maxres <- min(maxres, .Machine$integer.max)

    tune <- sqrt(qchisq(tune, p))

    if(n <= p + 1)
        stop(if (n <= p) "n <= p -- you can't be serious!" else "n == p+1  is too small sample size")
    if(n < 2 * p)
    { ## p+1 < n < 2p
        warning("n < 2 * p, i.e., possibly too small sample size")
    }

    ## Prepare for calling the Fortran function
    icent <- 1
    locat <- double(p)
    covmat <- matrix(0.0, p, p)
    storage.mode(covmat) <- "double"
    wk <- double(4*n+p)
    iwork <- integer(4*n+p)
    nresper <- 0
    w <- double(n)
    z <- double(n)

##  We are ignoring the 'center' parameter in the original function in 'robust'-
##      center = TRUE   - estimate center simultaneously with cov
##      center = FALSE  - assume the data are centered, i.e. center = rep(0,p)
##      center = vector of length p - sweep center from the data
##
##    if(length(center) == 1 && !center)
##        center <- rep(0, p)
##
##    if(length(center) > 1)
##    {
##        if(length(center) != p)
##            stop("the dimension of ", sQuote("center"), " does not match the ", "dimension of ", sQuote("x"))
##        x <- sweep(x, 2, center)
##        icent <- 0
##    }
##
    sdlist <- .Fortran("rlds",
                      n = as.integer(n),
                      p = as.integer(p),
                      nsamp = as.integer(nsamp),
                      x = as.double(x),
                      tune = as.double(tune),
                      wk = as.double(wk),
                      center = as.double(locat),
                      cov = covmat,
                      maxres = as.integer(maxres),
                      nresper = as.integer(nresper),
                      weights = as.double(w),
                      outlyingness = as.double(z),
                      icent = as.integer(icent),
                      iwork = as.integer(iwork),
                      PACKAGE = "rrcov")

## again skipping the predefined center
##
##    dist <- mahalanobis(x,
##        center = if(length(center) > 1) rep(0, p) else sdlist$center,
##        cov = sdlist$cov)

    center <- sdlist$center
    cov <- sdlist$cov
    mah <- mahalanobis(x, center=center, cov=cov)

    consistency.correction <- median(mah) / qchisq(.5, p)

##    cat("\nconsistency correction factor: ", consistency.correction, "\n")

    cov <- cov * consistency.correction
    mah <- mah / consistency.correction

##    if(length(center) > 1)
##        sdlist$center <- center

##    sdlist[c("cov", "center", "dist")]

    nms <- dimn[[2]]
    if(is.null(nms))
        nms <- paste("V", 1:p, sep = "")
    names(center) <- nms
    dimnames(cov) <- list(nms,nms)

    ans <- new("CovSde",
               call = xcall,
               iter=sdlist$nresper,                # actual subsamples performed
               crit=0,
               cov=cov,
               center=center,
               mah=mah,
               n.obs=n,
               wt=sdlist$weights,
               X = as.matrix(x),
               method="Stahel-Donoho estimator")
    ans
}
