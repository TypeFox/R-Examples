CovNAMcd <- function(x,
                   alpha=1/2,
                   nsamp=500,
                   seed=NULL,
                   trace=FALSE,
                   use.correction = TRUE,
                   impMeth = c("norm" , "seq", "rseq"),
                   control)
{

    ## Analize and validate the input parameters ...

    ## if a control object was supplied, take the option parameters from it,
    ## but if single parameters were passed (not defaults) they will override the
    ## control object.

    if(!missing(control)){
        defcontrol <- CovControlMcd()       # default control
        if(alpha == defcontrol@alpha)       alpha <- control@alpha
        if(nsamp == defcontrol@nsamp)       nsamp <- control@nsamp
        if(is.null(seed) || seed == defcontrol@seed)         seed <- control@seed
        if(trace == defcontrol@trace)       trace <- control@trace
        if(use.correction == defcontrol@use.correction) use.correction <- control@use.correction
    }

    impMeth <- match.arg(impMeth)

    ## prepare the call to covMcd() which will return an S3 object

    ## handle the case of nsamp="best" or "exact"
    iter <- if(is.numeric(nsamp)) nsamp else 0

    xcall <- match.call()
    mcd <- .covNAMcd(x=x, alpha=alpha, nsamp=nsamp, seed=seed, trace=trace, use.correction=use.correction, imputation=TRUE, impMeth=impMeth)
    ans <- new("CovNAMcd",
               call = xcall,
               iter=iter,
               crit=mcd$crit,
               cov=mcd$cov,
               center=mcd$center,
               n.obs=mcd$n.obs,
               X = mcd$X,
               method=mcd$method,
               best=mcd$best,
               alpha=mcd$alpha,
               quan=mcd$quan,
               raw.cov = mcd$raw.cov,
               raw.center = mcd$raw.center,
               raw.mah = mcd$raw.mah,
               raw.wt = mcd$raw.weights,
               raw.cnp2 = mcd$raw.cnp2,
               cnp2 = mcd$cnp2,
               singularity = mcd$singularity)
    ans
}
