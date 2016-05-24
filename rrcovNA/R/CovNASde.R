CovNASde <- function(x,
                   nsamp,
                   maxres,
                   tune=0.95,
                   eps=0.5,
                   prob=0.99,
                   impMeth = c("norm" , "seq", "rseq"),
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

    call <- match.call()


    ## drop all rows which contain only missings
    na.x <- rowSums(ifelse(is.na(x),1,0)) == ncol(x)
    ok <- !na.x
    x <- x[ok, , drop = FALSE]

    dimn <- dimnames(x)
    dx <- dim(x)
    n <- dx[1]
    p <- dx[2]
    if(p < 2)
        stop("Need at least 2 columns ")

##    s <- prelim.norm(x)                     # do preliminary manipulations
##    thetahat <- em.norm(s, showits=FALSE)   # find the mle estimates
##    rngseed(1234567)                        # set random number generator seed
##    ximp <- imp.norm(s, thetahat, x)        # impute missing data under the MLE
##    xx<-imp.norm(s, thetahat, x)

    ximp <- .imputation(x, impMeth = impMeth)
    sde <- CovSde(ximp, nsamp=nsamp, maxres=maxres, tune=tune, eps=eps, prob=prob, control=control)

    method="Stahel-Donoho estimator for incomplete data"
    ans <- new("CovNASde",
               call = call,
               iter=sde@iter,
               crit=0,
               cov=sde@cov,
               center=sde@center,
               n.obs=sde@n.obs,
               X = x,
               method=method)

    ans
}
