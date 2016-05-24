################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Calibration tests for count data based on proper scoring rules
### Reference: Wei and Held (2014), Test, 23, 787-805
###
### Copyright (C) 2015 Sebastian Meyer
### $Revision: 1424 $
### $Date: 2015-07-15 12:14:54 +0200 (Mit, 15. Jul 2015) $
################################################################################


## perform a calibration test given observations x
## with Poisson (size = NULL) or NegBin predictions

calibrationTest.default <- function (x, mu, size = NULL,
                                     which = c("dss", "logs", "rps"),
                                     tolerance = 1e-4, method = 2, ...)
{
    stopifnot(x >= 0, mu > 0, is.null(size) || size > 0)
    
    ## calculate scores
    which <- match.arg(which)
    score <- do.call(which, args = alist(x = x, mu = mu, size = size))

    ## calculate z-statistic
    z <- calibrationZ(score, mu, size, which, tolerance, method)

    ## calculate two-sided p-value
    p <- 2 * pnorm(-abs(z))
    
    ## construct an object of class "htest"
    res <- list(
        method = paste0("Calibration Test for Count Data (based on ",
            toupper(which), ")"),
        data.name = deparse(substitute(x)),
        statistic = c("z" = z),
        parameter = c("n" = length(x)),
        p.value = p
    )
    class(res) <- "htest"
    res
}

## compute the calibration z-statistic given the computed scores
calibrationZ <- function (score, mu, size = NULL,
                          which = c("dss", "logs", "rps"),
                          tolerance = 1e-4, method = 2)
{
    stopifnot(method %in% 1:2)
    ## expectation and variance of score for given predictive distribution
    EV <- score_EV(mu, size, tolerance, which)
    ## calculate the z-statistic
    z <- do.call(paste0("zScore", method),
                 args = alist(score, EV[[1L]], EV[[2L]]))
    z
}

## compute the calibration z-statistic and p-value
## from a set of scores and their null expectations and variances
zScore1 <- function (score, E0, V0)
{
    n <- length(score)

    ## emean <- mean(E0)
    ## varmean <- sum(V0) / n^2
    ## (mean(score) - emean) / sqrt(varmean)
    
    sum(score - E0) / sqrt(sum(V0))
}

## alternative z-statistic Z*
zScore2 <- function (score, E0, V0)
{
    n <- length(score)
    sum((score - E0) / sqrt(V0)) / sqrt(n)
}
