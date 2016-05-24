# Constructor function for a bds identity matrix
#  The first arg will become the dimnames
#
bdsI <- function(id, blocksize) {
    n <- length(id)
    if (n==1 && is.integer(id) && id >0) {
        # like diag(), we allow a simple count
        bdsmatrix(blocksize=rep(1,id), blocks=rep(1., id))
        }
    else {
        if (missing(blocksize)) {
            bdsmatrix(blocksize=rep(1,n), blocks=rep(1., n),
                      dimnames=list(id,id))
            }
        else {
            if (sum(blocksize) != length(id)) stop("Inconsitent arguments")
            temp <- sum(blocksize*(blocksize+1)/2)
            x <- bdsmatrix(blocksize=blocksize, blocks=rep(0., temp),
                           dimnames=list(id,id))
            diag(x) <- rep(1.0, length(id))
            x
            }
        }
    }
