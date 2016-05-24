CovMcd <- function(x,
                   raw.only = FALSE,
                   alpha=control@alpha,
                   nsamp = control@nsamp,
                   scalefn=control@scalefn, maxcsteps=control@maxcsteps,
                   initHsets=NULL, save.hsets=FALSE,
                   seed=control@seed,
                   trace=control@trace,
                   use.correction=control@use.correction,
                   control=CovControlMcd(), ...)
{
    ## MM 2014-11-04: add '...' so we are compatible to covMcd() extensions in robustbase
    ## Analize and validate the input parameters ...
    ## if a control object was supplied, take the option parameters from it,
    ## but if single parameters were passed (not defaults) they will override the
    ## control object.
    ## prepare the call to covMcd() which will return an S3 object

    ## handle the case of nsamp="best" or "exact"
    iter <- if(is.numeric(nsamp)) nsamp else 0

    xcall <- match.call()
    mcd <- covMcd(x=x, raw.only=raw.only, alpha=alpha, nsamp=nsamp,
                  scalefn=control@scalefn,
                  maxcsteps=control@maxcsteps,
                  initHsets = NULL, save.hsets = FALSE,
                  seed=seed, trace=trace, use.correction=use.correction, ...)

    new("CovMcd",
        call= xcall,
        iter=iter,
        crit=mcd$crit,
        cov=mcd$cov,
        center=mcd$center,
        n.obs=mcd$n.obs,
        mah = mcd$mah,
        wt = mcd$mcd.wt,
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
}
