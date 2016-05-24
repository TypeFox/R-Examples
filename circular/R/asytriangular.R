#############################################################
#                                                           #
#   dasytriangular function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: July, 20, 2011                                    #
#   Copyright (C) 2011 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################
#Mardia pag. 52
dasytriangular <- function(x, rho) {
   if (rho < 0 | rho > 1/pi)
      stop("'rho' must be between 0 and 1/pi")
   x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
   attr(x, "class") <- attr(x, "circularp") <-  NULL

   DasytriangularRad(x, rho)
}

DasytriangularRad <- function(x, rho) {
   d <- (1 + rho * (x - pi))/(2 * pi)
   return(d)
}

