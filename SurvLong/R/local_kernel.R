#******************************************************************************#
# Kernel functions                                                             #
#******************************************************************************#
#                                                                              #
# Inputs                                                                       #
#                                                                              #
#  t              an object of class numeric.                                  #
#                 the time point(s) at which the kernel is to be calculated    #
#                                                                              #
#  h              an object of class  numeric.                                 #
#                 the kernel bandwidth                                         #
#                                                                              #
#  kType          an object of class character indicating the type of          #
#                 smoothing kernel to use in the estimating equation.          #
#                 Must be one of \{"epan", "uniform", "gauss"\}, where         #
#                 "epan" is the Epanechnikov kernel and "gauss" is the         #
#                 Gaussian kernel.                                             #
#                                                                              #
# Outputs                                                                      #
#                                                                              #
#  An object of class numeric.                                                 #
#                                                                              #
#******************************************************************************#
local_kernel <- function(t, h, kType){

  if( kType == "epan" ) {

    res <- epanechnikov(t/h)/h

  } else if( kType == "uniform" ) {

    res <- uniform(t/h)/h

  } else if( kType == "gauss" ) {

    res <- gauss(t/h)/h

  } else {

    stop("unsupported kernel")

  }

  return(res)
}

epanechnikov <- function(t){

  tst <- (-1.0 <= t) & (t <= 1.0 )

  kt <- 0.75*( 1.0 - t*t )
  kt[ !tst ] <- 0.0

  return(kt)
}

uniform <- function(t){

  tst <- (-1.0 <= t) & (t <= 1.0 )
  kt <- t
  kt[  tst ] <- 0.5
  kt[ !tst ] <- 0.0

  return(kt)
}

gauss <- function(t){

  kt <- exp(-t*t*0.5)/sqrt(2.0*pi)

  return(kt)

}
