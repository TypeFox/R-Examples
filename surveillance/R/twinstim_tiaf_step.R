################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Step function implementation for temporal interaction
###
### Copyright (C) 2014 Sebastian Meyer
### $Revision: 735 $
### $Date: 2014-01-31 22:52:57 +0100 (Fre, 31 Jan 2014) $
################################################################################


tiaf.step <- function (knots, maxRange = Inf, nTypes = 1, validpars = NULL)
{
    knots <- sort(unique(as.vector(knots,mode="numeric")))
    stopifnot(knots > 0, is.finite(knots), isScalar(maxRange), maxRange > knots)
    nknots <- length(knots)             # = number of parameters (per type)
    knotsmax <- c(knots, maxRange)
    nknotsmax <- nknots + 1L
    allknots <- c(0, knots, maxRange)
    nallknots <- length(allknots)
    allknotsInf <- unique(c(allknots, Inf)) # ensure Inf as last element

    stopifnot(isScalar(nTypes <- as.integer(nTypes)), nTypes > 0L)
    npars <- nknots * nTypes
    .parintwidths <- rep.int(c(diff.default(knotsmax)), nTypes)
    .parintwidths[is.infinite(.parintwidths)] <- -1
    ##<- in case maxRange=Inf, last interval width will always be multiplied by
    ##   0 and should give 0, but Inf would produce NaN, so we just set it to -1

    ## the step function is right-continuous, intervals are [a,b)
    g <- if (nTypes > 1) {
        heights <- function (logvals) { # get matrix of type-specific heights
            dim(logvals) <- c(nknots, nTypes)
            rbind(1, exp(logvals), 0, deparse.level=0)
        }
        function (t, logvals, types)
            heights(logvals)[(types-1)*nallknots + .bincode(t, allknotsInf, right=FALSE)]
    } else {
        function (t, logvals, types = NULL)
            c(1,exp(logvals),0)[.bincode(t, allknotsInf, right=FALSE)]
    }

    G <- if (nTypes > 1) {
        typeheights <- function (logvals, type) # vector of type-specific heights
            c(1, exp(logvals[(type-1)*nknots+seq_len(nknots)]))
        as.function(c(alist(t=, logvals=, types=), substitute({
            mapply(function (t, type) {
                knots2t <- c(0, pmin.int(knots, t), TMAX)
                sum(typeheights(logvals, type) * diff.default(knots2t))
            }, t, types, SIMPLIFY=TRUE, USE.NAMES=FALSE)
        }, list(TMAX = if (is.finite(maxRange))
                quote(min(t,maxRange)) else quote(t)))))
    } else {
        ## function (t, logvals, types = NULL) {
        ##     vapply(t, function (t) {
        ##         knots2t <- c(0, pmin.int(knots, t), min(t, maxRange))
        ##         sum(c(1,exp(logvals)) * diff.default(knots2t))
        ##     }, 0, USE.NAMES=FALSE)      # vapply is faster than sapply
        ## }
        as.function(c(alist(t=, logvals=, types = NULL), substitute({
            ##omtk <- outer(t, knots, pmin.int), bare-bone implementation:
            omtk <- pmin.int(rep.int(knots, rep.int(L <- length(t), nknots)), t)
            dim(omtk) <- c(L, nknots)
            .colSums(apply(cbind(0, omtk, TMAX, deparse.level=0),
                           1L, diff.default) * c(1,exp(logvals)),
                     nknotsmax, L)
        }, list(TMAX = if (is.finite(maxRange))
                quote(pmin.int(t,maxRange)) else quote(t)))))
    }

    ## the derivative is simply the height corresponding to (t, type) and is 0
    ## outside this interval/type
    deriv <- function (t, logvals, types) {
        whichvals <- .bincode(t, knotsmax, right=FALSE)
        fixedheight <- is.na(whichvals)
        ##<- intervals number 1 and 'nallknots' don't correspond to parameters
        whichvals <- whichvals + (types-1)*nknots # select type parameter
        whichvals[fixedheight] <- 0
        ## we do a bare-bone implementation of relevant parts of
        ## deriv <- outer(whichvals, seq_len(npars), "==") * rep(exp(logvals), each=L)
        repL <- rep.int(L <- length(t), npars)
        Y <- rep.int(seq_len(npars), repL) # column index
        Z <- rep.int(exp(logvals), repL)   # value
        ##<- 6x faster than rep(..., each=L)
        res <- (Y == whichvals) * Z
        dim(res) <- c(L, npars)
        res
    }
    ## only tiny modification necessary for nTypes == 1
    if (nTypes == 1) {
        body(deriv)[[grep("types", body(deriv))]] <- NULL
        formals(deriv)["types"] <- list(NULL)
    }

    Deriv <- deriv
    body(Deriv) <- as.call(append(as.list(body(Deriv)), expression(
        partwidth <- t - knots[whichvals]
        ), after=2L))
    body(Deriv)[[grep("whichvals[fixedheight]", body(Deriv), fixed=TRUE)]] <-
        quote(whichvals[fixedheight] <- partwidth[fixedheight] <- 0)
    body(Deriv) <- as.call(append(as.list(body(Deriv)), expression(
        W <- rep.int(.parintwidths, repL)
        ), after=grep("Z <-", body(Deriv))))
    body(Deriv)[[grep("res <-", body(Deriv))]] <- if (nTypes == 1) {
        quote(res <- ((Y < whichvals | t >= maxRange) * W + (Y == whichvals) * partwidth) * Z)
    } else {
        quote(res <- ((Y > (types-1)*nknots & (Y < whichvals | t >= maxRange)) * W +
                      (Y == whichvals) * partwidth) * Z)
    }

    ## Done
    res <- list(g = g, G = G,
                deriv = deriv, Deriv = Deriv,
                ## FIXME: simulate = simulate,
                npars = npars, validpars = validpars)
    attr(res, "knots") <- knots
    attr(res, "maxRange") <- maxRange
    res
}
