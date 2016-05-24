###
### $Id: midwt.R 22 2014-06-20 20:59:33Z plroebuck $
### Computes the inverse discrete wavelet transform x for a 1D or
###           2D input signal y using the scaling filter h.
###

##-----------------------------------------------------------------------------
midwt <- function(y, h, L) {
    .Call("do_midwt",
          as.matrix(y),
          as.vector(h),
          as.integer(L),
          PACKAGE="rwt")
}

