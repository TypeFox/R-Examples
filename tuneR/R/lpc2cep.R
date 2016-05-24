# This code is based on the Matlab implementations of PLP and Rasta
# feature calculations by Daniel P. W. Ellis of Columbia University /
# International Computer Science Institute.  For more details, see:
# http://www.ee.columbia.edu/~dpwe/resources/matlab/rastamat/

lpc2cep <- function(a, nout=nrow(a)){

    if(!(is.numeric(a) && is.matrix(a)))
      stop("'a' has to be a numeric matrix")

    if(!(nout==as.integer(nout) && nout > 0))
        stop("'modelorder' has to be a positive integer")

    arow <- nrow(a)
    acol <- ncol(a)

    order <- arow - 1

    mc <- matrix(0, nrow=nout, ncol=acol)

    # Code copied from HSigP.c: LPC2Cepstrum

    # First cep is log(Error) from Durbin
    mc[1,] = -log(a[1,])

    # Renormalize lpc A coeffs
    a <- a / matrix(rep(a[1,], arow), nrow=arow, byrow=TRUE)

    for(n in 2:nout){
        s <- 0
        for(m in 2:n){
            s <- s + (n-m) * a[m,] * mc[n-m+1,]
        }
        mc[n,] <- -(a[n,] + s/(n-1))
    }
    features <- mc

    return(features)
}

