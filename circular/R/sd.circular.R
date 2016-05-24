
sd <- function(x, ...) UseMethod("sd")

sd.default <- function(x, na.rm = FALSE, ...) stats::sd(x=x, na.rm=na.rm)

sd.data.frame <- function(x, ...) {
    sapply(x, sd, ...)
}

##############################################################
#                                                            #
#   sd.circular function                                     #
#   Author: Claudio Agostinelli and Jean-Olivier Irisson     #
#   Email: claudio@unive.it                                  #
#   Date: June, 24, 2011                                     #
#   Copyright (C) 2011 Claudio Agostinelli                   #
#                                                            #
#   As defined in                                            #
#   Mardia, KV. Statistics of directional data. 1972         #
#   Formula actually taken from                              #
#   Zar, JH. Biostatistical analysis. 2010, 26.5, p 617      #
#                                                            #
#   Version 0.3                                              #
##############################################################

sd.circular <- function (x, na.rm=FALSE, ...)  {
  if (is.matrix(x)) {
	  # NB: matrices cannot be handled by a method because a matrix of circular data would have "circular" as its first class
    apply(x, 2, sd.circular, na.rm=na.rm)

  } else {
    # Remove missing values
    if (na.rm) {
      x <- x[!is.na(x)]
    }

    # Checks
    if (length(x) == 0) {
      warning("No observations (at least after removing missing values)")
      return(NA)
    }

    # Possibly set and then get the circular attributes of the input data
    if (!is.circular(x)) {
      x <- circular(x)
    }
    datacircularp <- circularp(x)

    # Compute the standard deviation
    x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
    attr(x, "class") <- attr(x, "circularp") <- NULL
    s <- SdCircularRad(x=x)

    # Determine circular attributes of interest for the result,
    # based on the attributes of the data and on control.circular
##    dc <- control.circular
##    if (is.null(dc$type))
##      dc$type <- datacircularp$type
##    if (is.null(dc$units))
##     dc$units <- datacircularp$units

##    # Convert the standard deviation into the appropriate circular class
##    s <- conversion.circular(circular(s), units=dc$units, type=dc$type, template="none", modulo="2pi", zero=0, rotation="counter")
    return(s)
  }
}

SdCircularRad <- function(x) {
  rbar <- RhoCircularRad(x)
  circsd <- sqrt(-2*log(rbar))
  return(circsd)
}
