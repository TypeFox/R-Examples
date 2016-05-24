# This code is based on the Matlab implementations of PLP and Rasta
# feature calculations by Daniel P. W. Ellis of Columbia University /
# International Computer Science Institute.  For more details, see:
# http://www.ee.columbia.edu/~dpwe/resources/matlab/rastamat/

deltas <- function(x, w=9){

    if(!(is.numeric(x) && is.matrix(x)))
      stop("'x' has to be a numeric matrix")

    if(!(w==as.integer(w) && w > 0))
      stop("'w' has to be a positive integer")

    nr <- nrow(x)
    nc <- ncol(x)

    # Define windows shape
    hlen <- floor(w/2)
    w <- 2 * hlen + 1
    win <- seq(hlen, -hlen, -1)

    # Pad data
    xx <- matrix(c(rep(x[,1], hlen), x, rep(x[,nc], hlen)), nrow=nr)

    # Delta filtering alog rows
    d <-  t(apply(xx, 1, function(x) {convolve(x, win, conj=FALSE, type="open")[-(1:(length(win)-1))]}))

    # Trim edges
    d <- d[,2 * hlen + 1:nc,drop=FALSE]

    return(d)
}

