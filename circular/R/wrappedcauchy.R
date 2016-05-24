
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   rwrappedcauchy function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: August, 10, 2006                                  #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-3                                           #
#############################################################

rwrappedcauchy <- function(n, mu = circular(0), rho = exp(-1), control.circular=list()) {
   if (is.circular(mu)) {
      datacircularp <- circularp(mu)
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
   
   mu <- conversion.circular(mu, units="radians", zero=0, rotation="counter")
   attr(mu, "class") <- attr(mu, "circularp") <-  NULL  
   result <- RwrappedcauchyRad(n, mu, rho)
   result <- conversion.circular(circular(result), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)

   return(result)
}

RwrappedcauchyRad <- function(n, mu, rho) {
    if (rho == 0)
    result <- runif(n, 0, 2 * pi)
    else if (rho == 1)
         result <- rep(mu, n)
    else {
       scale <-  - log(rho)
       result <- rcauchy(n, mu, scale) %% (2 * pi)
    }
    return(result)
}

#############################################################
#                                                           #
#   dwrappedcauchy function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: May, 22, 2006                                     #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2                                             #
#############################################################

dwrappedcauchy <- function(x, mu=circular(0), rho=exp(-1)) {
    if (rho < 0 | rho > 1)
        stop("rho must be between 0 and 1")
    x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
    mu <- conversion.circular(mu, units="radians", zero=0, rotation="counter")

    attr(x, "class") <- attr(x, "circularp") <-  NULL
    attr(mu, "class") <- attr(mu, "circularp") <-  NULL    

    DwrappedcauchyRad(x, mu, rho)
}

DwrappedcauchyRad <- function(x, mu, rho) {
    d <- (1 - rho^2)/((2 * pi) * (1 + rho^2 - 2 * rho * cos(x - mu)))
    return(d)
}
