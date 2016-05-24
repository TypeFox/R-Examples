###
### $Id: mdwt.R 22 2014-06-20 20:59:33Z plroebuck $
### Computes the discrete wavelet transform y for a 1D or 2D input
###          signal x using the scaling filter h.
###

##-----------------------------------------------------------------------------
mdwt <- function(x, h, L) {
    .Call("do_mdwt",
          as.matrix(x),
          as.vector(h),
          as.integer(L),
          PACKAGE="rwt")
}

