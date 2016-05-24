################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Student (t) kernel f(s) = (||s||^2+sigma^2)^-d
### This is a reparametrization of the t-kernel; For d=1, this is the kernel of
### the Cauchy density with scale sigma; in Geostatistics, a correlation
### function of this kind is known as the Cauchy model.
###
### Copyright (C) 2013-2014 Sebastian Meyer
### $Revision: 711 $
### $Date: 2014-01-27 21:51:24 +0100 (Mon, 27 Jan 2014) $
################################################################################


siaf.student <- function (nTypes = 1, validpars = NULL)
{
    nTypes <- as.integer(nTypes)
    stopifnot(length(nTypes) == 1L, nTypes > 0L)

    ## for the moment we don't make this type-specific
    if (nTypes != 1) stop("type-specific shapes are not yet implemented")

    ## helper expression, note: logpars=c(logscale=logsigma, logd=logd)
    tmp <- expression(
        logsigma <- logpars[[1L]],  # used "[[" to drop names
        logd <- logpars[[2L]],
        sigma <- exp(logsigma),
        d <- exp(logd)
        )

    ## spatial kernel
    f <- function (s, logpars, types = NULL) {}
    body(f) <- as.call(c(as.name("{"),
        tmp,
        expression(s2 <- .rowSums(s^2, nrow(s), 2L)),
        expression((s2+sigma^2)^-d)
    ))

    ## numerically integrate f over a polygonal domain
    F <- function (polydomain, f, logpars, type = NULL, ...)
        .polyCub.iso(polydomain$bdry, intrfr.student, logpars, #type,
                     center=c(0,0), control=list(...))
    
    ## fast integration of f over a circular domain
    ## is not relevant for this heavy-tail kernel since we don't use
    ## 'effRange', and usually eps.s=Inf
    ##Fcircle <- function (r, logpars, type = NULL) {}

    ## derivative of f wrt logpars
    deriv <- f
    body(deriv)[[length(body(deriv))]] <- # assignment for return value of f
        substitute(fvals <- x, list(x=body(deriv)[[length(body(deriv))]]))
    body(deriv) <- as.call(c(as.list(body(deriv)), expression(
        derivlogsigma <- -2*d*sigma^2 * fvals / (s2+sigma^2),
        derivlogd <- log(fvals) * fvals,
        cbind(derivlogsigma, derivlogd, deparse.level = 0)
        )))

    ## Numerical integration of 'deriv' over a polygonal domain
    Deriv <- function (polydomain, deriv, logpars, type = NULL, ...)
    {
        res.logsigma <- .polyCub.iso(polydomain$bdry,
                                     intrfr.student.dlogsigma, logpars, #type,
                                     center=c(0,0), control=list(...))
        res.logd <- .polyCub.iso(polydomain$bdry,
                                 intrfr.student.dlogd, logpars, #type,
                                 center=c(0,0), control=list(...))
        c(res.logsigma, res.logd)
    }

    ## simulation from the kernel
    simulate <- siaf.simulatePC(intrfr.student)

    ## set function environments to the global environment
    environment(f) <-  environment(deriv) <- .GlobalEnv
    ## in F, Deriv, and simulate we need access to the intrfr-functions
    environment(F) <- environment(Deriv) <- environment(simulate) <-
        getNamespace("surveillance")

    ## return the kernel specification
    list(f=f, F=F, deriv=deriv, Deriv=Deriv, simulate=simulate,
         npars=2L, validpars=validpars)
}


## integrate x*f(x) from 0 to R (vectorized)
intrfr.student <- function (R, logpars, types = NULL)
{
    sigma <- exp(logpars[[1L]])
    d <- exp(logpars[[2L]])
    if (d == 1) {
        log(R^2+sigma^2) / 2 - log(sigma)
    } else {
        ( (R^2+sigma^2)^(-d+1) - (sigma^2)^(-d+1) ) / (2-2*d)
    }
}

## integrate x * (df(x)/dlogsigma) from 0 to R (vectorized)
intrfr.student.dlogsigma <- function (R, logpars, types = NULL)
{
    sigma <- exp(logpars[[1L]])
    d <- exp(logpars[[2L]])
    sigma^2 * ( (R^2+sigma^2)^-d - sigma^(-2*d) )
}

## integrate x * (df(x)/dlogd) from 0 to R (vectorized)
intrfr.student.dlogd <- function (R, logpars, types = NULL)
{
    sigma <- exp(logpars[[1L]])
    d <- exp(logpars[[2L]])
    if (d == 1) {
        log(sigma)^2 - log(R^2+sigma^2)^2 / 4
    } else { # thanks to Maple 17
        primitive <- function (x) {
            x2ps2 <- x^2 + sigma^2
            (d*(d-1)*log(x2ps2) + d) / (2*(d-1)^2 * (x2ps2)^(d-1))
        }
        primitive(R) - primitive(0)
    }
}
