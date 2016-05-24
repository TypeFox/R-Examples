################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Expectation and variance of proper scoring rules for Poisson and NegBin
### Reference: Wei and Held (2014), Test, 23, 787-805
###
### Copyright (C) 2013-2014 Wei Wei, 2015 Sebastian Meyer
### $Revision: 1512 $
### $Date: 2015-11-05 10:11:03 +0100 (Don, 05. Nov 2015) $
################################################################################

## wrapper function calling the necessary "EV" function for the selected score
score_EV <- function (mu, size = NULL, tolerance = 1e-4,
                      which = c("dss", "logs", "rps"))
{
    which <- match.arg(which)
    if (which == "dss")
        return(dss_EV(mu, size))

    ## for "logs" and "rps", the EV function only works with a single prediction
    ## -> apply to each mu (size)
    res <- if (is.null(size)) { # Poisson
        vapply(X = mu,
               FUN = paste0(which, "_EV_1P"),
               FUN.VALUE = c(E = 0, V = 0),
               tolerance = tolerance,
               USE.NAMES = FALSE)
    } else { # NegBin
        mapply(FUN = paste0(which, "_EV_1NB"),
               mu = mu, size = size,
               MoreArgs = list(tolerance = tolerance),
               SIMPLIFY = TRUE, USE.NAMES = FALSE)
    }
    ## 'res' has dimension 2 x length(mu)
    list(E = res[1L,], V = res[2L,])
}


##########################
### Dawid-Sebastiani Score
##########################

dss_EV <- function (mu, size = NULL)
{
    sigma2 <- if (is.null(size)) mu else mu * (1 + mu/size)
    E <- 1 + log(sigma2)
    V <- if (is.null(size)) {
        2 + 1/sigma2
    } else {
        2 + 6/size + 1/sigma2
    }
    list(E = E, V = V)
}


#####################
### Logarithmic Score
#####################

## for a single Poisson prediction
logs_EV_1P <- function (mu, tolerance = 1e-4) # tolerance is in absolute value
{
    ## use the same kmax for expectation and variance -> shared computations
    ## K2 is always a bit larger than K1, so we use K2
    kmax <- if (mu^3 < tolerance/.Machine$double.eps/2) {
        ## we can calculate K2 from Theorem 1 (b)
        qpois(1 - tolerance/(mu^3 + 6*mu^2 + 7*mu + 1), lambda = mu) + 3
    } else { # very high quantile (e.g., 1 - 1e-16) would yield Inf
        mu + 10 * sqrt(mu)
    }
    kseq <- seq_len(kmax)

    ## compute values required by both E and V
    fseq <- dpois(kseq, lambda = mu)
    logfactseq <- lfactorial(kseq)
    
    ## expectation
    E <- if (mu > tolerance^(-1/4)) { # fast version for "large" mu
        ## approximation error is of order 1/mu^4
        0.5 + 0.5*log(2*pi*mu) - 1/12/mu - 1/24/mu^2 - 19/360/mu^3
    } else {
        ##kmax1 <- qpois(1 - tolerance/(mu^2 + 3*mu + 1), lambda = mu) + 2
        seqq1 <- fseq * logfactseq
        mu * (1-log(mu)) + sum(seqq1)
    }
    
    ## variance (does it converge to 0.5 as mu -> Inf ?)
    seqq2 <- (logfactseq - kseq * log(mu))^2 * fseq
    V <- sum(seqq2) - (E - mu)^2
    
    c(E = E, V = V)
}

## for a single NegBin prediction
logs_EV_1NB <- function (mu, size, tolerance = 1e-4)
{
    ## TODO: replace simple kmax by formulae from the paper
    kmax <- qnbinom(1-tolerance/10, mu = mu, size = size) + 5
    kseq <- 0:kmax

    ## compute values required by both E and V
    fseq <- dnbinom(kseq, mu = mu, size = size)
    lgammaseq <- lbeta(kseq + 1L, size) + log(kseq + size)
    
    ## expectation
    seqq1 <- lgammaseq * fseq
    E <- sum(seqq1) - size*log(size) - mu*log(mu) + (mu+size)*log(mu+size)

    ## variance
    con2 <- E - size * log(1 + mu/size)
    seqq2 <- (lgammaseq + kseq * log(1 + size/mu))^2 * fseq
    V <- sum(seqq2) - con2^2
    ## check against formulation in the paper (Equation 11):
    ## con2paper <- E + size*log(size) - size*log(size+mu) - lgamma(size)
    ## seqq2paper <- (-lgamma(kseq+size) + lgamma(kseq+1L) + kseq*log(1+size/mu))^2 * fseq
    ## Vpaper <- sum(seqq2paper) - con2paper^2
    ## => V and Vpaper are only identical for kmax -> Inf
    
    c(E = E, V = V)
}


############################
### Ranked Probability Score
############################

## for a single Poisson prediction
rps_EV_1P <- function (mu, tolerance = 1e-4) # tolerance is in absolute value
{
    ## expectation
    if (requireNamespace("gsl", quietly = TRUE)) {
        ## faster and more accurate implementation (works for larger mu)
        E <- mu * gsl::bessel_I0_scaled(2*mu, give=FALSE, strict=TRUE) +
            mu * gsl::bessel_I1_scaled(2*mu, give=FALSE, strict=TRUE)
    } else {
        E <- mu * besselI(2*mu, 0, expon.scaled = TRUE) +
            mu * besselI(2*mu, 1, expon.scaled = TRUE)
        if (identical(E, 0)) {
            ## R's besselI() works fine for mu <= 50000 (on my .Machine)
            ## but returns 0 (buffer overflow) for larger arguments
            warning("'mu' is too large for besselI(), install package \"gsl\"")
            return(c(E = NA_real_, V = NA_real_))
        }
    }

    ## variance
    kmax <- max(qpois(1 - tolerance/(10*mu^2 + mu), lambda = mu) + 2,
                8)  # cf. Theorem 2 (a)
    kseq <- 0:kmax
    fseq <- dpois(kseq, lambda = mu)
    Fseq <- cumsum(fseq)  # = ppois(kseq, lambda = mu)
    psiseq <- (kseq - mu) * (2*Fseq - 1) + 2*mu * fseq
    seqq <- psiseq^2 * fseq
    V <- sum(seqq) - 4 * E^2
    
    c(E = E, V = V)
}

## for a single NegBin prediction
rps_EV_1NB <- function (mu, size, tolerance = 1e-4)
{
    ## determine kmax for Var0(RPS), which is always > kmax for E0(RPS),
    ## cf. Theorem 2 (c), here corrected (1-) and simplified
    l5 <- (mu + 1)^2 + 1
    kmax2 <- max(qnbinom(1-tolerance/l5, mu = mu*(1+2/size), size = size+2) + 2,
                 8)  
    ## the other listed terms seem to be always smaller than the first one:
    ## qnbinom(1-tolerance/l5, mu = mu, size = size)
    ## qnbinom(1-tolerance/l5, mu = mu*(1+1/size), size = size+1) + 1
    kseq2 <- 0:kmax2
    fseq2 <- dnbinom(kseq2, mu = mu, size = size)
    Fseq2 <- cumsum(fseq2)  # = pnbinom(kseq2, mu = mu, size = size)

    ## expectation
    ghgz_part <- mu * (1 + mu/size)
    ghgz <- 4 * ghgz_part / size
    E <- if (ghgz < 1 && requireNamespace("gsl", quietly = TRUE)) {
        ghgz_part * gsl::hyperg_2F1(1+size, 0.5, 2, -ghgz, give = FALSE, strict = TRUE)
    } else {
        kmax1 <- max(qnbinom(1-tolerance/mu, mu = mu*(1+1/size), size = size+1) + 1,
                     8)  # cf. Theorem 2 (b)
        kseq1 <- seq_len(kmax1)
        seqq1 <- vapply(
            X = kseq1, # we could use kmax2 (> kmax1) also here
            FUN = function (i) fseq2[i+1L] * sum((i:1) * fseq2[seq_len(i)]),
            FUN.VALUE = 0, USE.NAMES = FALSE)
        sum(seqq1)
    }
    
    ## variance
    psiseq <- kseq2 * (2 * Fseq2 - 1) +
        mu * (1 - 2 * pnbinom(kseq2 - 1, mu = mu + mu/size, size = size + 1))
    seqq <- psiseq^2 * fseq2
    V <- sum(seqq) - 4 * E^2
    
    c(E = E, V = V)
}
