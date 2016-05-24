##
##  o d r l i n r e g . R  Orthogonal Distance Regression
##


# Linear orthogonal distance regression method
odregress <- function(x, y) {
    stopifnot(is.numeric(x), is.numeric(y))
    
    Z <- cbind(x, y)
    n <- nrow(Z)      # no. of data points
    m <- ncol(Z) - 1  # no. of independent variables

    meanZ <- repmat(apply(Z, 2, mean), n, 1)
    svdZ <- svd(Z - meanZ)
    V <- svdZ$v

    a <- -V[1:m, m+1] / V[m+1, m+1]
    b <- mean(Z %*% V[, m+1]) / V[m+1, m+1]

    # Fitted values
    yfit <- cbind(x, 1) %*% c(a, b)
    resd <- y - yfit

    # orthogonal distance
    normal <- V[, m+1]
    err <- abs((Z - meanZ) %*% normal)
    ssq <- sum(err^2)

    return( list(coeff = c(a, b), ssq = ssq, err = err, 
                 fitted = yfit, resid = resd, normal = normal) )
}
