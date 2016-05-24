
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   rho.circular function                                   #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: May, 26, 2006                                     #
#   Version: 0.3-1                                          #
#                                                           #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#############################################################

rho.circular <- function(x, na.rm=FALSE) {
   if (na.rm) 
       x <- x[!is.na(x)]
   if (any(is.na(x))) {
       warning("No observations (at least after removing missing values)")
       return(NA)
   }
   x <- conversion.circular(x, units="radians")
   RhoCircularRad(x)
}

RhoCircularRad <- function(x) {
   n <- length(x)   
   sinr <- sum(sin(x))
   cosr <- sum(cos(x))
   result <- sqrt(sinr^2 + cosr^2)/n
   return(result)
}
