
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   trigonometric.moment function                           #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: August, 10, 2006                                  #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.3-3                                           #
#############################################################

trigonometric.moment <- function(x, p = 1, center = FALSE, control.circular=list()) {
   x <- unlist(x)
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

   # Handling missing values
   x <- na.omit(x)
   if ((n <- length(x))==0) {
      warning("No observations (at least after removing missing values)")
      return(NULL)
   }      
   x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
   attr(x, "class") <- attr(x, "circularp") <-  NULL
    
   res <- TrigonometricMomentRad(x, p, center)
   mu.p <- conversion.circular(circular(res[1]), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
   result <- list(mu=mu.p, rho=res[2], cos=res[3], sin=res[4], p=res[5], n=res[6], call=match.call())
   return(result)
}

TrigonometricMomentRad <- function(x, p, center) {
    center <- as.numeric(center)
    n <- length(x)
    sinr <- sum(sin(x))
    cosr <- sum(cos(x))
    circmean <- atan2(sinr, cosr)
    sin.p <- sum(sin(p * (x - circmean * center)))/n
    cos.p <- sum(cos(p * (x - circmean * center)))/n
    mu.p <- atan2(sin.p, cos.p)
    rho.p <- sqrt(sin.p^2 + cos.p^2)
    result <- c(mu.p, rho.p, cos.p, sin.p, p, n)
    return(result)
}

