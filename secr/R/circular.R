################################################################################## package 'secr'
## circular.R
## "circular error probable"
## last changed 2011 06 12; 2013-04-19; 2013-04-24; 2013-05-11
################################################################################

circular.r <- function (p = 0.95, detectfn = 0, sigma = 1, detectpar = NULL, hazard = TRUE, ...) {

    ## translate character detectfn to numeric code
    if (is.character(detectfn))
        detectfn <- detectionfunctionnumber(detectfn)

    if (detectfn %in% c(0,2,3,14,16)) {
        ## if input is a named list
        if (!is.null(detectpar)) {
            sigma <- detectpar$sigma
        }
        detectpar <- list(sigma = sigma)
    }
    else
        if (is.null(detectpar))
            stop ("require detectpar, in list format, except for ",
                  "halfnormal, exponential and uniform")
    detectpar$g0 <- 1  ## always
    detectpar$lambda0 <- 1  ## 2013-04-19
    truncate <- detectpar$truncate
    if (is.null(truncate))
         truncate <- Inf
    OK <- truncate == Inf
    detectpar <- detectpar[parnames(detectfn)]  ## correct order
    pars <- unlist(detectpar)
    cutval <- ifelse (detectfn %in% c(9,10,11), detectpar$cutval, NA)
    scale <- spatialscale (detectpar, detectfn) ## see pdot.R; assumes cutval in detectpar

    ## use formula for halfnormal
    if (OK & (((detectfn == 0) & !hazard) | ((detectfn == 14) & hazard))) {
        (-2*log(1-p))^0.5 * sigma
    }
    else if (OK & (((detectfn == 2) & !hazard) | ((detectfn == 16) & hazard))) {
        fnr <- function (r, this.p) {
            1 - (r/sigma + 1) * exp(-r/sigma) - this.p
        }
        getroot <- function (p) uniroot(fnr, c(0,200*scale), this.p = p)$root
        sapply(p, getroot)
    }
    ## uniform is dead easy
    else if (OK & ((detectfn == 3) & !hazard)) {
        p^0.5 * sigma
    }
    ## otherwise integrate
    else {
        dfn <- getdfn (detectfn)
        rdfn <- function (r, pars, cutval)  {
            haz <- dfn(r, pars, cutval)
            if (hazard) haz <- -log(1-haz)
            haz[!is.finite(haz)] <- 0
            r * haz
        }
        I1 <- integrate (rdfn, 0, truncate, pars, cutval, ...)$value

        fnr <- function (r, this.p) {
            I2 <- integrate (rdfn, 0, min(r,truncate), pars, cutval, ...)$value
            I2 / I1 - this.p
        }
        getroot <- function (p) uniroot(fnr, c(0,200*scale), this.p = p)$root
        sapply(p, getroot)
    }
}


circular.p <- function (r = 1, detectfn = 0, sigma = 1, detectpar = NULL, hazard = TRUE, ...) {

    ## convert character detectfn to numeric code
    if (is.character(detectfn))
        detectfn <- detectionfunctionnumber(detectfn)

    if (detectfn %in% c(0,2,3,14,16)) {
        ## if input is a named list
        if (!is.null(detectpar)) {
            sigma <- detectpar$sigma
        }
        detectpar <- list(sigma = sigma)
    }
    else
        if (is.null(detectpar))
            stop ("require detectpar, in list format, except for ",
                  "halfnormal, exponential and uniform")
    detectpar$g0 <- 1       ## always
    detectpar$lambda0 <- 1  ## always
    truncate <- detectpar$truncate
    if (is.null(truncate))
         truncate <- Inf
    OK <- truncate == Inf
    detectpar <- detectpar[parnames(detectfn)]  ## correct order
    pars <- unlist(detectpar)
    cutval <- ifelse (detectfn %in% c(9,10,11), detectpar$cutval, NA)

    scale <- spatialscale (detectpar, detectfn) ## see pdot.R; assumes cutval in detectpar

    ## use formula for halfnormal
    if (OK & (((detectfn == 0) & !hazard) | ((detectfn == 14) & hazard))) {
        1 - exp(-(r/sigma)^2 / 2)
    }
    else if (OK & (((detectfn == 2) & !hazard) | ((detectfn == 16) & hazard))) {
        1 - (r/sigma + 1) * exp(-r/sigma)
    }
    ## uniform is dead easy
    else if (OK & (detectfn == 3) & !hazard) {
        (r/sigma)^2
    }
    ## otherwise integrate
    else {
        dfn <- getdfn(detectfn)
        rdfn <- function (r, pars, cutval) {
            haz <- dfn(r, pars, cutval)
            if (hazard) haz <- -log(1-haz)
            haz[!is.finite(haz)] <- 0
            r * haz
        }
        I1 <- integrate (rdfn, 0, truncate, pars, cutval, ...)$value

        fnr <- function (r) {
            I2 <- integrate (rdfn, 0, min(r,truncate), pars, cutval, ...)$value
            I2 / I1
        }
        sapply(r, fnr)
    }
}

# plot(seq(0,5,0.1), circular.p(seq(0,5,0.1), detectfn=0), type='l', xlab='radius', ylab='p')
# lines(seq(0,5,0.1),circular.p(seq(0,5,0.1), detectfn=1, detectpar=list(sigma=1, z=4)), col='blue')
# lines(seq(0,5,0.1),circular.p(seq(0,5,0.1), detectfn=2), col='red')



