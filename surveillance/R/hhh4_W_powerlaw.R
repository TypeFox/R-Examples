################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Parametric power-law specification for neighbourhood weights in hhh4()
###
### Copyright (C) 2012-2015 Sebastian Meyer
### $Revision: 1353 $
### $Date: 2015-06-03 21:46:45 +0200 (Mit, 03. Jun 2015) $
################################################################################


### Construct weight matrix wji according to the Zeta-distribution with respect
### to the orders of neighbourhood (in nbmat, as e.g. obtained from nbOrder()),
### optionally fulfilling rowSums(wji) = 1
## As a formula (for j != i, otherwise wji = 0):
## wji = pzeta(oji; d, maxlag) / sum_k pzeta(ojk; d, maxlag)
## Here, oji = oij is the order of nb of i and j,
## and pzeta(o; d, m) = o^-d / sum_{r=1}^m r^-d is the Zeta-distribution
## on 1:m (also called Zipf's law).
## Special cases: maxlag >= max(nbmat) yields the weights
## wji = oji^-d / sum_k ojk^-d
## and maxlag=1 yields the classical weights wji=1/nj.

zetaweights <- function (nbmat, d = 1, maxlag = max(nbmat), normalize = FALSE)
{
    ## raw (non-normalized) zeta-distribution on 1:maxlag
    zeta <- c(0, seq_len(maxlag)^-d)  # first 0 is for lag 0 (i.e., diag(nbmat))

    ## replace order by zetaweight of that order
    wji <- zeta[nbmat + 1L]           # results in vector
    wji[is.na(wji)] <- 0              # for lags > maxlag

    ## set dim and names
    dim(wji) <- dim(nbmat)
    dimnames(wji) <- dimnames(nbmat)

    if (normalize) normalizeW(wji) else wji
}



### powerlaw weights
## in the non-truncated case, i.e. maxlag = max(nbmat),
## the raw powerlaw weights are defined as w_ji = o_ji^-d,
## and with (row-)normalization we have    w_ji = o_ji^-d / sum_k o_jk^-d

W_powerlaw <- function (maxlag, normalize = TRUE, log = FALSE,
                        initial = if (log) 0 else 1)
{
    if (missing(maxlag)) {
        stop("'maxlag' must be specified (e.g. maximum neighbourhood order)")
        ## specifying 'maxlag' in zetaweights is actually optional since it has
        ## the default value max(nbmat). however, repeatedly asking for this
        ## maximum would be really inefficient.
    } else stopifnot(isScalar(maxlag))

    ## main function which returns the weight matrix
    weights.call <- call("zetaweights",
                         quote(nbmat), quote(d), maxlag, normalize)
    weights <- as.function(c(alist(d=, nbmat=, ...=), call("{", weights.call)),
                           envir=.GlobalEnv)
    if (log) { # the parameter d is interpreted on log-scale
        ## we prepend the necessary conversion d <- exp(d)
        body(weights) <- as.call(append(as.list(body(weights)),
                                        quote(d <- exp(d)), after=1))
    }
    
    ## construct derivatives with respect to "d" (or log(d), respectively)
    dweights <- d2weights <- as.function(c(alist(d=, nbmat=, ...=), quote({})),
                                         envir=.GlobalEnv)
    weights.call[[5L]] <- FALSE         # normalize separately
    header <- c(
        if (log) quote(d <- exp(d)),    # such that d is again on original scale
        substitute(Wraw <- weights.call, list(weights.call=weights.call)),
        if (normalize) expression(
            nUnits <- nrow(Wraw),
            norm <- .rowSums(Wraw, nUnits, nUnits)
            ),
        expression(
            # Wraw == 0 means o = 0 (diagonal) or o > maxlag  =>  deriv = 0
            is.na(Wraw) <- Wraw == 0,  # set to NA since we will take the log
            logo <- -log(Wraw)/d       # = log(nbmat) with NA's at Wraw == 0 
            ),
        if (normalize) quote(W <- Wraw / norm) else quote(W <- Wraw)
        )
    footer <- expression(deriv[is.na(deriv)] <- 0, deriv)

    ## first derivative
    tmp1 <- expression(
        ## in surveillance < 1.9-0, 'norm' and 'tmpnorm' were based on 'nbmat',
        ## which is incorrect for the truncated case maxlag < max(nbmat)
        tmpnorm <- .rowSums(Wraw * -logo, nUnits, nUnits, na.rm=TRUE) / norm,
        tmp1 <- logo + tmpnorm
        )
    deriv1 <- if (normalize) {
        expression(deriv <- W * -tmp1)
    } else expression(deriv <- W * -logo)
    body(dweights) <- as.call(c(as.name("{"),
            header,
            if (normalize) tmp1,
            deriv1,
            if (log) expression(deriv <- deriv * d), # this is the non-log d
            footer
        ))

    ## second derivative
    body(d2weights) <- as.call(c(as.name("{"),
            header,
            if (normalize) {
                c(tmp1, expression(
                    tmp2 <- .rowSums(Wraw * logo^2, nUnits, nUnits,
                                     na.rm=TRUE) / norm - tmpnorm^2,
                    deriv <- W * (tmp1^2 - tmp2)
                    ))
            } else expression(deriv <- W * logo^2),
            if (log) c(
                do.call("substitute", list(deriv1[[1L]], list(deriv=as.name("deriv1")))),
                expression(deriv <- deriv * d^2 + deriv1 * d) # this is the non-log d
                ),
            footer
        ))

    ## return list of functions
    list(w=weights, dw=dweights, d2w=d2weights, initial=initial)
}
