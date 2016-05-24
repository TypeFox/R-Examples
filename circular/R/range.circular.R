#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################


#############################################################
#                                                           #
#   range.circular function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: May, 06, 2011                                     #
#   Copyright (C) 2011 Claudio Agostinelli                  #
#                                                           #
#   Version 0.5                                             #
#############################################################

range.circular <- function(x, test = FALSE, na.rm=FALSE, finite=FALSE, control.circular=list(), ...) {
  if (finite) 
    x <- x[is.finite(x)]
  if (na.rm) 
    x <- x[!is.na(x)]
  else {
    if (any(is.na(x))) {
       x <- circular(NA)
       return(x)
    }
  }
  if (is.circular(x)) {
    datacircularp <- circularp(x)     
  } else {
    datacircularp <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")
  }

  dc <- control.circular
  if (is.null(dc$type))
    dc$type <- datacircularp$type
  if (is.null(dc$units))
    dc$units <- datacircularp$units
  if (is.null(dc$template))
    dc$template <- datacircularp$template
  if (is.null(dc$modulo))
    dc$modulo <- datacircularp$modulo
  if (is.null(dc$zero))
    dc$zero <- datacircularp$zero
  if (is.null(dc$rotation))
    dc$rotation <- datacircularp$rotation
   
  x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
  attr(x, "class") <- attr(x, "circularp") <- NULL

  result <- RangeCircularRad(x, test)
   
  if (test) {
    result$range <- conversion.circular(x=circular(result$range, template=dc$template, rotation='counter'), units=dc$units, type=dc$type, modulo="asis", zero=NULL)
  } else {
    result <- conversion.circular(x=circular(result, template=dc$template, rotation='counter'), units=dc$units, type=dc$type, modulo="asis", zero=NULL)
  }
  return(result)
}

RangeCircularRad <- function(x, test=TRUE) {
   x <- sort(x %% (2*pi))
   n <- length(x)
   spacings <- c(diff(x), x[1] - x[n] + 2*pi)
   range <- 2*pi - max(spacings)
   if (test == TRUE) {
       stop <- floor(1/(1 - range/(2*pi)))
       index <- c(1:stop)
       sequence <- ((-1)^(index - 1)) * exp(lgamma(n + 1) - lgamma(index + 1) - lgamma(n - index + 1)) * (1 - index * (1 - range/(2 * pi)))^(n - 1)
       p.value <- sum(sequence)
       result <- list(range=range, p.value=p.value)
   } else {
       result <- range
   }
   return(result)
}
