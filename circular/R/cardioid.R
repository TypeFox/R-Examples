
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   rcardioid function                                      #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: August, 10, 2006                                  #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-3                                           #
#############################################################

rcardioid <- function(n, mu=circular(0), rho=0, control.circular=list()) {
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
   res <- RcardioidRad(n, mu, rho)
   res <- conversion.circular(circular(res), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
   return(res)
}

RcardioidRad <- function(n, mu, rho) {
   if (rho < -0.5 | rho > 0.5)
       stop("rho must be between -0.5 and 0.5")        
   i <- 1
   result <- rep(0, n)
   while (i <= n) {
      x <- runif(1, 0, 2 * pi)
      y <- runif(1, 0, (1 + 2 * rho)/(2 * pi))
      f <- (1 + 2 * rho * cos(x - mu))/(2 * pi)
      if (y <= f) {
         result[i] <- x
         i <- i + 1
      }
   }
   return(result)
}

#############################################################
#                                                           #
#   dcardioid function                                      #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: May, 31, 2006                                     #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.3-1                                           #
#############################################################

dcardioid <- function(x, mu=circular(0), rho=0) {
    x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
    mu <- conversion.circular(mu, units="radians", zero=0, rotation="counter")

    attr(x, "class") <- attr(x, "circularp") <-  NULL
    attr(mu, "class") <- attr(mu, "circularp") <-  NULL    
  
    DcardioidRad(x, mu, rho)
}

DcardioidRad <- function(x, mu=0, rho=0) {
    if (rho < -0.5 | rho > 0.5)
        stop("rho must be between -0.5 and 0.5")
    d <- (1 + 2 * rho * cos(x - mu))/(2 * pi)
    return(d)
}
