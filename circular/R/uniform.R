#############################################################
#                                                           #
#   rcircularuniform function                               #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: November, 06, 2013                                #
#   Copyright (C) 2013 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1-1                                           #
#############################################################

rcircularuniform <- function(n, control.circular=list()) {
   datacircularp <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")
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
   
   un <- RuniformRad(n)
   un <- conversion.circular(circular(un), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
   return(un)
}

RuniformRad <- function(n) {
   un <- stats::runif(n=n, min=0, max=2*pi)
   return(un)
}

#############################################################
#                                                           #
#   dcircularuniform function                               #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: March, 31, 2009                                   #
#   Copyright (C) 2009 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

dcircularuniform <- function (x) {
  x <- rep(1/(2*pi), length(x))
  return(x)
}

