# This code is based on the Matlab implementations of PLP and Rasta
# feature calculations by Daniel P. W. Ellis of Columbia University /
# International Computer Science Institute.  For more details, see:
# http://www.ee.columbia.edu/~dpwe/resources/matlab/rastamat/

dolpc <- function(x, modelorder=8){

    if(!(is.numeric(x) && is.matrix(x)))
      stop("'x' has to be a numeric matrix")

    if(!(modelorder == as.integer(modelorder) && modelorder > 0))
        stop("'modelorder' has to be a positive integer")

    nbands <- nrow(x)

    # px <- planFFT(2*nbands-2)
    # Calculate autocorrelation 
    r <- apply(rbind(x, x[seq(nbands-1, 2, -1),]), 2,
            function(y) Re(fft(y, inverse=TRUE))/length(y))
            # function(y) Re(IFFT(y, plan=px)))
    # First half only
    r <- r[1:nbands,,drop=FALSE]

    # Find LPC coeffs by Levinson-Durbin
    levcoef <- levinson(x=r, p=modelorder)

    # Normalize each poly by gain
    y <- t(levcoef$a) / matrix(rep(levcoef$v, modelorder+1), nrow=modelorder+1,
            byrow=TRUE)

    return(y)
}

