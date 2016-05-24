#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[invGamma.R] by DSB Fre 04/07/2008 11:36 (CEST) on daniel@puc.home>
##
## Description:
## Functions for the inverse gamma distribution.
##
## History:
## 04/07/2008   copy from thesis function collection.
#####################################################################################

dinvGamma <- function (x, a, b, log = FALSE, normalize = TRUE)
{
    ret <- - (a + 1) * log (x) - b / x
    if (normalize)
        ret <- ret + a * log (b) - lgamma (a)
    if (log)
        return (ret)
    else
        return (exp (ret))
}

##

pinvGamma <- function (q, a, b, lower.tail = TRUE, log.p = FALSE)
{
    pgamma (q = 1 / q, shape = a, rate = b, lower.tail = !lower.tail, log.p = log.p)
}

##

qinvGamma <- function (p, a, b, lower.tail = TRUE, log.p = FALSE)
{
    1 / qgamma (p = p, shape = a, rate = b, lower.tail = !lower.tail, log.p = log.p)
}

##

rinvGamma <- function (n, a, b)
{
    1 / rgamma (n, shape = a, rate = b)
}

