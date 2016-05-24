#############################################################
#                                                           #
#   quantile.circular function                              #
#   Author: Claudio Agostinelli and Alessandro Gagliardi    #
#   Email: claudio@unive.it                                 #
#   Date: August, 26, 2013                                  #
#   Copyright (C) 2013 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-1                                           #
#############################################################


quantile.circular <- function(x, probs = seq(0, 1, 0.25), na.rm=FALSE, names = TRUE, type = 7, ...) {
   if (na.rm)
      x <- x[!is.na(x)]
   if (length(x)==0) {
      warning("No observations (at least after removing missing values)")
      return(NULL)
   }
   if(probs < 0 || probs > 1){
      warning("'probs' outside [0,1]")
      return(NULL)   
   }
   if (is.circular(x)) {
      dc <- circularp(x)
   } else {
      dc <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")
   }
   x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
   attr(x, "class") <- attr(x, "circularp") <-  NULL
   circquantile <- QuantileCircularRad(x=x, probs=probs, names=names, type=type, ...)
   circquantile <- conversion.circular(circular(drop(circquantile)), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)  
   return(circquantile)
}


## quantile.default becomes quantile
## 20130826
QuantileCircularRad <- function(x, probs = seq(0, 1, 0.25), names = TRUE, type = 7, ...) {
   circularmedian <- MedianCircularRad(x)
   if(is.na(circularmedian))
      return(rep(NA,length(probs)))
   attr(circularmedian, "medians") <- NULL
   tx <- (x-circularmedian)%%(2*pi)
   tx <- MinusPiPlusPiRad(tx)
   circularQuantile <- quantile(x=tx, probs=probs, names=names, type=type, ...)
   return((circularQuantile + circularmedian)%%(2*pi))   
}
