
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   rtriangular function                                    #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: May, 29, 2006                                     #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-2                                           #
#############################################################

rtriangular <- function(n, rho, control.circular=list()) {
    dc <- control.circular
    theta <- RtriangularRad(n, rho)
    theta <- conversion.circular(circular(theta), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
    return(theta)
}

RtriangularRad <- function(n, rho) {
   if (rho < 0 | rho > 4/pi^2)
      stop("'rho' must be between 0 and 4/pi^2")
   u <- matrix(c(runif(n)), ncol = 1)
   get.theta <- function(u, rho) {
      if (u < 0.5) {
         a <- pi * rho
         b <-  - (4 + pi^2 * rho)
         c <- 8 * pi * u
         theta1 <- ( - b + sqrt(b^2 - 4 * a * c))/(2 * a)
         theta2 <- ( - b - sqrt(b^2 - 4 * a * c))/(2 * a)
         min(theta1, theta2)
      } else {
         a <- pi * rho
         b <- 4 - 3 * pi^2 * rho
         c <- (2 * pi^3 * rho) - (8 * pi * u)
         theta1 <- ( - b + sqrt(b^2 - 4 * a * c))/(2 * a)
         theta2 <- ( - b - sqrt(b^2 - 4 * a * c))/(2 * a)
         max(theta1, theta2)
      }
   }
   theta <- apply(u, 1, get.theta, rho)
   theta[theta > pi] <- theta[theta > pi] - 2 * pi
   return(theta)
}

#############################################################
#                                                           #
#   dtriangular function                                    #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: May, 24, 2006                                     #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2                                             #
#############################################################

dtriangular <- function(x, rho) {
   if (rho < 0 | rho > 4/pi^2)
      stop("'rho' must be between 0 and 4/pi^2")
   x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
   attr(x, "class") <- attr(x, "circularp") <-  NULL

   DtriangularRad(x, rho)
}

DtriangularRad <- function(x, rho) {
   d <- (4 - pi^2 * rho + 2 * pi * rho * abs(pi - x))/(8 * pi)
   return(d)
}

