##
##  h a m p e l . R  MAD Outlier in Time Series
##


hampel <- function (x, k, t0 = 3)
{
    #   x:  vector or time series
    #   k:  window [x_(i-k),...,x_i,...,x_(i+k)]
    n   <- length(x)
    y   <- x         # corrected x vector
    ind <- c()       # indices of outliers

    L  <- 1.4826     # constants for normal distributions
    # t0 <- 3        # Pearson's 3 sigma edit rule

    # we don't look at outliers at the end parts of x !
    for ( i in (k+1):(n-k) ) {
        x0 <- median( x[(i-k):(i+k)] )
        S0 <- L * median( abs(x[(i-k):(i+k)] - x0) )
        if ( abs(x[i]-x0) > t0 * S0 ) {
            y[i] <- x0
            ind  <- c(ind, i)
        }
    }
    # return a list with 2 components
    list(y=y, ind=ind)
}
