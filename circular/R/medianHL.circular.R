#############################################################
#                                                           
#   median.circular function                                  
#   Author: Claudio Agostinelli and Alessandro Gagliardi
#   E-mail: claudio@unive.it                                
#   Date: October, 12, 2012                                  
#   Version: 0.4                                          
#                                                           
#   Copyright (C) 2012 Claudio Agostinelli and Alessandro Gagliardi
#                                                           
#############################################################

medianHL.circular <- function(x, na.rm=FALSE, method=c("HL1","HL2","HL3"), prop=NULL) {
   method <- match.arg(method)
   if (!is.null(prop))
      if (prop <= 0 | prop >=1)
         stop("'prop' is outside (0,1)")
   if (na.rm)
      x <- x[!is.na(x)]
   if (length(x)==0) {
      warning("No observations (at least after removing missing values)")
      return(NULL)
   }
   if (is.circular(x)) {
      dc <- circularp(x)
   } else {
      dc <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")
   }
   x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
   attr(x, "class") <- attr(x, "circularp") <-  NULL
   circmedian <- MedianHLCircularRad(x, method, prop)
   circmedian <- conversion.circular(circular(drop(circmedian)), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
   attr(circmedian, "medians") <- conversion.circular(circular(drop(MinusPiPlusPiRad(attr(circmedian, "medians")))), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)  
   attr(attr(circmedian, "medians"), "class") <- attr(attr(circmedian, "medians"), "circularp") <-  NULL
   return(circmedian)
}

MedianHLCircularRad <- function(x, method, prop) {
   x <- x%%(2*pi)
   x <- MinusPiPlusPiRad(x)
   n <- length(x)
   mediancirc = NA
	methods <- c("HL2","HL1","HL3")
   if (is.null(prop))
	{
   	mediancirc <- .C("MedianHLCircularRad",x=as.double(x),y=as.double(x),n=as.integer(n),whichMethod=as.integer(which(methods==method) - 1),result=as.double(0))$result
   }
	else
	{
		mediancirc <- .C("MedianHLCircularPropRad",x=as.double(x),n=as.integer(n),whichMethod=as.integer(which(methods==method) - 1),prop=as.double(prop),result=as.double(0))$result
   }
   return(mediancirc)
}
