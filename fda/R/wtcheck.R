wtcheck = function(n, wtvec=NULL) {
#  WTVEC is either a vector or a matrix, and WTCHECK applies
#  a number of tests to confirm that will served as a weight
#  vector or matrix for smoothing data.  n is the required length
#  or order of WTVEC.  The functions returns WTVEC, ONEWT indicating
#  whether WTVEC contains only one's, and MATWT indicating that
#  WTVEC is in fact a square positive definite matrix.

#  Last modified 9 July 2011 by Jim Ramsay

#  check n

if (n != round(n)) stop("n is not an integer.")
if (n < 1)         stop("n is less than 1.")

#  check wtvec

if (!is.null(wtvec)) {
    dimw = dim(as.matrix(wtvec))
    if (any(is.na(as.vector(wtvec)))) stop("WTVEC has NA values.")
    if (all(dimw == n)) {
        #  WTVEC is a matrix of order n
        onewt = FALSE
        matwt = TRUE
        #  check weight matrix for being positive definite
        wteig = eigen(wtvec)$values
        if (any(is.complex(wteig))) stop("Weight matrix has complex eigenvalues.")
        if (min(wteig) <= 0)        stop("Weight matrix is not positive definite.")
    } else {
        #  WTVEC is treated as a vector
        if ((length(dimw) > 1 && dimw[1] > 1 && dimw[2] > 1) || length(dimw) > 2) {
            stop ("WTVEC is neither a vector nor a matrix of order n.")
        }
        wtvec = as.matrix(wtvec)
        if (length(wtvec) == 1) {
            wtvec = wtvec*matrix(1,n,1)
        }
        if (length(wtvec) != n) {
            stop("WTVEC of wrong length")
        }
        if (min(wtvec) <= 0) stop("Values in WTVEC are not positive.")
        onewt = FALSE	
        matwt = FALSE
    }
} else {
    wtvec = matrix(1,n,1)
    onewt = TRUE
    matwt = FALSE
}

return(list(wtvec=wtvec, onewt=onewt, matwt=matwt))

}

