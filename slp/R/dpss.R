#
#     Originally written by Karim Rahim based on Percival and Walden (1993), and
#     updated to use LAPACK. Makes use of technique found in David Thomson's F77 code for
#     reducing the tridiagonal matrix in half, based on Slepian technical memo,
#     Bell Labs (1977). 
#
#     Changes made by Wesley Burr, 2013. Eigenvalues not included in formulation,
#     as they are not needed, and speed is a priority.
#
#     If you wish to report bugs please contact the maintainer:
# 
#     Wesley Burr
#     <wesley.burr@gmail.com>

##########################################################
#
#  .dpss
# 
#  Generates k orthogonal Discrete Prolate Spheroidal
#  Sequences (dpss) using the tridiagonal method. See
#  Slepian (1978), page 1379 and Percival and Walden (1993),
#  Chapter 8.4. 
#
#  Very similar to 'dpss' in pkg:multitaper, but without 
#  eigenvalues.
#
##########################################################

.dpss <- function(n, k, nw) {

    stopifnot(n >= 1, nw/n >0, nw/n < 0.5, k >= 1)
    # if(k > 2 * nw) { warning("dpss: Generating sequences with _very_ poor in-band properties. K > 2NW.") }

    ## If k is passed in as floating point, the cast to 
    ## as.integer() in the Fortran call does not quite work properly
    ## - wburr, 2013
    if(!is.integer(k)) {
      k<-as.integer(floor(k));
    } 

    ## 'eigen' is of appropriate length for LAPACK functions.
    ## This call uses LAPACK functions in place of the
    ## EISPACK functions referenced in Percival and Walden (1993). 
    ## Code ported to F90 by krahim. 
    out <- .Fortran("dpss", as.integer(n), as.integer(k),
              as.double(nw), v=double(n*k), eigen=double(k),
              PACKAGE='slp')
    out <- matrix(data=out$v, nrow=n, ncol=k, byrow=FALSE)

    class(out) <- "dpss"
    return(out)
}

