################################################################
#This function computes the density of a multivariate          #
#normal distribution. I could have used dmvnorm, but           #
#to avoid spending time in checking whether the covaiance      #
#matrix is positive definite, I got the part that only computes#
#the density, since when I use this function at the IMIS       # #algorithm, the covariance matrix is surely positive definite  #
################################################################

dmnorm_mod <- function (x, mean = rep(0, d), d, varcov, invvarcov) 
    {
    x <- matrix(x, 1, d)
    n <- 1
    X <- t(matrix(x, nrow = n, ncol = d)) - mean
    Q <- apply((invvarcov %*% X) * X, 2, sum)
    logDet <- sum(logb(abs(diag(qr(varcov)[[1]]))))
    logPDF <- as.vector(Q + d * logb(2 * pi) + logDet)/(-2)
    return(exp(logPDF))
    }