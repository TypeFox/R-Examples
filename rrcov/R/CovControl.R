CovControlMcd <- function (alpha=0.5,
                           nsamp=500,
                           scalefn=NULL,
                           maxcsteps=200,
                           seed=NULL,
                           trace=FALSE,
                           use.correction=TRUE)
{
    new("CovControlMcd", alpha = alpha,
                         nsamp = nsamp,
                         scalefn=scalefn,
                         maxcsteps=maxcsteps,
                         seed = seed,
                         trace = trace,
                         use.correction = use.correction)
}

setMethod("restimate", "CovControlMcd", function(obj, x, ...)
    CovMcd(x, control = obj, ...)
)

CovControlMest <- function (r = 0.45,
                            arp = 0.05,
                            eps=1e-3,
                            maxiter=120
                           )
{
    new("CovControlMest", r = r, arp = arp, eps=eps, maxiter=maxiter)
}

setMethod("restimate", "CovControlMest", function(obj, x, ...)
    CovMest(x, control = obj, ...)
)

CovControlOgk <- function (niter=2,
                           beta=0.90,
                           mrob=NULL,
                           vrob=.vrobGK,
                           smrob="scaleTau2",
                           svrob="gk"
                           )
{
    new("CovControlOgk", niter=niter, beta=beta, mrob=mrob, vrob=vrob, smrob=smrob, svrob=svrob)
}

setMethod("restimate", "CovControlOgk", function(obj, x, ...)
    CovOgk(x, control = obj, ...)
)

CovControlMve <- function (alpha=0.5,
                           nsamp=500,
                           seed=NULL,
                           trace=FALSE)
{
    new("CovControlMve", alpha = alpha,
                         nsamp = nsamp,
                         seed = seed,
                         trace = trace)
}

setMethod("restimate", "CovControlMve", function(obj, x, ...)
    CovMve(x, control = obj, ...)
)


## Moved to AllClassess - should precede the definition of class CovMMest,
## sinse used in the prototype
##
##CovControlSest <- function (bdp=0.5,
##                            arp=0.1,
##                            eps=1e-5,
##                            maxiter=120,
##                            nsamp=500,
##                            seed=NULL,
##                            trace=FALSE,
##                            tolSolve=1e-14,
##                            method="sfast")
##{
##    new("CovControlSest", bdp=bdp, arp=arp, eps=eps, maxiter=maxiter,
##        nsamp=nsamp, seed=seed, trace=trace, tolSolve=tolSolve, method=method)
##}
##

setMethod("restimate", "CovControlSest", function(obj, x, ...)
    CovSest(x, control = obj, ...)
)

CovControlSde <- function  (nsamp=0,
                            maxres=0,
                            tune=0.95,
                            eps=0.5,
                            prob=0.99,
                            seed=NULL,
                            trace=FALSE,
                            tolSolve=1e-14)
{
    new("CovControlSde", nsamp=nsamp, maxres=maxres, tune=tune, eps=eps, prob=prob,
        seed=seed, trace=trace, tolSolve=tolSolve)
}

setMethod("restimate", "CovControlSde", function(obj, x, ...)
    CovSde(x, control = obj, ...)
)

CovControlMMest <- function (bdp=0.5,
                            eff=0.95,
                            maxiter=50,
                            sest=CovControlSest(),
                            trace=FALSE,
                            tolSolve=1e-7)
{
    new("CovControlMMest", bdp=bdp, eff=eff, maxiter=maxiter, sest=sest,
        trace=trace, tolSolve=tolSolve)
}

setMethod("restimate", "CovControlMMest", function(obj, x, ...)
    CovMMest(x, control = obj, ...)
)
