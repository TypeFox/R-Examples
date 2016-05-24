# 
# Matrix multiplication for symmetric block diagonal (bds) matrices
#
bdsmult <- function(x, y) {
    dy <- dim(y)
    dx <- dim(x)
    ldy <- length(dy)
    if (ldy!=2) dy <- c(length(y), 1) # y is a vector
    if (dx[2] != dy[1]) 
	stop("Number of columns of x should be the same as number of rows of y")

    # Do the multiplication in C code.  Y is replaced by the result
    #  (Since x is a square matrix, the result is the same size as y)
    nblock <- length(x@blocksize)
    temp <- .C("bdsmatrix_prod", 
	       as.integer(nblock),
	       as.integer(x@blocksize),
	       as.integer(dy),
	       as.double(x@blocks),
	       as.double(x@rmat),
	       as.double(x@offdiag),
	       temp = double(dy[1]),
	       itemp= integer(max(1,x@blocksize)),
	       y =   as.double(y))
    z <- matrix(temp$y, nrow=dx[1])

    # Create dimnames for the result, using the dimnames of the input args
    dnx <- dimnames(x)
    dny <- dimnames(y)
    if(!is.null(dnx) || !is.null(dny)) {
	dnz <- list(NULL, NULL)
	if(!is.null(dnx))
	    dnz[1] <- dnx[1]
	if(!is.null(dny))
	    dnz[2] <- dny[2]
	dimnames(z) <- dnz
        }
    z
    }

setMethod("%*%", signature(x='bdsmatrix', y='matrix'), bdsmult)
setMethod("%*%", signature(x='bdsmatrix', y='numeric'), bdsmult)

#
# This allows for multiplication in the other direction
#
setMethod("%*%", signature(x='matrix', y='bdsmatrix'),
    function(x, y) {
	t(y%*% t(x))
	})

setMethod("%*%", signature(x='numeric', y='bdsmatrix'),
    function(x, y) {
	t(y%*% x)
	})

