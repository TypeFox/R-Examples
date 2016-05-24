#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[rmvt.R] by DSB Fre 02/10/2009 10:49 (CEST)>
##
## Description:
## Helper function for sampling from multivariate normal distribution.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 02/10/2009   don't rely on rmvnorm package
#####################################################################################

rmvt <- function (n, sigma = diag (2), mu = rep (0, 2), df = 1)
{
    ## simulate multivariate normal, using the symmetric root of sigma
    ev <- eigen(sigma, symmetric = TRUE)   
    if (! all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
        warning("sigma is numerically not positive definite")
    }
    
    ret <- ev$vectors %*% diag(sqrt(ev$values), length(ev$values)) %*% t(ev$vectors) 
    ret <- matrix(rnorm(n * ncol(sigma)), nrow = n) %*% ret

    ## then go to multivariate Student
    ret <- ret / sqrt(rchisq(n, df) / df)

    ## make the location shift
    ret <- sweep (ret, 2, mu, "+")

    ## return the samples matrix
    ret
}
