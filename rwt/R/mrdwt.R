###
### $Id: mrdwt.R 22 2014-06-20 20:59:33Z plroebuck $
### Computes the redundant discrete wavelet transform y
###           for a 1D or 2D input signal.
###

##-----------------------------------------------------------------------------
mrdwt <- function(x, h, L) {
    .Call("do_mrdwt",
          as.matrix(x),
          as.vector(h),
          as.integer(L),
          PACKAGE="rwt")
}

