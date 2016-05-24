#############################################################
#                                                           #
#   angular.deviation function                              #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: June, 24, 2011                                    #
#   Copyright (C) 2011 Claudio Agostinelli                  #
#                                                           #
#   Version 0.3                                             #
#############################################################

angular.deviation <- function (x, na.rm=FALSE)  {
  if (is.matrix(x)) {
    apply(x, 2, angular.deviation, na.rm=na.rm)
  } else {
    if (na.rm) 
      x <- x[!is.na(x)]
    x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
    attr(x, "class") <- attr(x, "circularp") <-  NULL
    AngularDeviationRad(x=x)
  }
}

AngularDeviationRad <- function(x) {
  rbar <- RhoCircularRad(x)
  circvar <- sqrt(2*(1-rbar))
  return(circvar)
}


