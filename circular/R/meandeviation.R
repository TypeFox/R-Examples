#############################################################
#                                                           #
#   meandeviation function                                  #
#   Author: Claudio Agostinelli and Alessandro Gagliardi    #
#   Email: claudio@unive.it                                 #
#   Date: August, 29, 2012                                  #
#   Copyright (C) 2012 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1                                             #
#############################################################

meandeviation <- function (x, na.rm=FALSE)  {
  if (is.matrix(x)) {
    apply(x, 2, meandeviation, na.rm=na.rm)
  } else {
    if (na.rm) 
      x <- x[!is.na(x)]
    x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
    attr(x, "class") <- attr(x, "circularp") <-  NULL
    MeanDeviationRad(x=x)
  }
}

MeanDeviationRad <- function(x) {
  circmedian <- MedianCircularRad(x)
##  circmedian <- attr(circmedian, 'medians')[1]
  meandev <- pi - mean(abs(pi-abs(MinusPiPlusPiRad(x-circmedian))))
  return(meandev)
}
