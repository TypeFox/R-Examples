
#############################################################
#                                                           #
#   daxialvonmises function                                 #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: May, 24, 2006                                     #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2                                             #
#############################################################

daxialvonmises <- function (x, mu, kappa, l=2) {
   if (l<=0)
      stop("'l' must be non negative")
   x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
   mu <- conversion.circular(mu, units="radians", zero=0, rotation="counter")
    
   attr(x, "class") <- attr(x, "circularp") <-  NULL
   attr(mu, "class") <- attr(mu, "circularp") <-  NULL

   DaxialvonmisesRad(x, mu, kappa, l)
}

DaxialvonmisesRad <- function(x, mu, kappa, l=2) {
    d <- l/(2 * pi * besselI(x = kappa, nu = 0, expon.scaled = TRUE)) * (exp(cos(l*(x - mu))-1))^kappa
    return(d)
}
