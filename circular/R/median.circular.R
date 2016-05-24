### This is necessary since stats::median do not have ... argument
### Work around suggested by Kurt but does not work.
###median <- function(x, na.rm, ...) UseMethod("median")
###median.default <- function(x, na.rm, ...) stats::median.default(x, na.rm)

#############################################################
#                                                           
#   median.circular function                                  
#   Author: Claudio Agostinelli and Alessandro Gagliardi
#   E-mail: claudio@unive.it                                
#   Date: September, 11, 2012                                  
#   Version: 0.4                                          
#                                                           
#   Copyright (C) 2012 Claudio Agostinelli and Alessandro Gagliardi
#                                                           
#############################################################

median.circular <- function(x, na.rm=FALSE) {
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
  x <- conversion.circular(x, units="radians")
  attr(x, "class") <- attr(x, "circularp") <-  NULL
  circmedian <- MedianCircularRad(x)
  circmedian <- conversion.circular(circular(drop(circmedian), template=dc$template, zero=dc$zero, rotation=dc$rotation), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
  attr(circmedian, "medians") <- conversion.circular(circular(drop(attr(circmedian, "medians")), template=dc$template, zero=dc$zero, rotation=dc$rotation), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)  
  attr(attr(circmedian, "medians"), "class") <- attr(attr(circmedian, "medians"), "circularp") <-  NULL
  return(circmedian)
}

MedianCircularRad <- function(x)
{
	n <- length(x)
	res <- .C("MedianCircularRad",x=as.double(x),n=as.integer(n),result=as.double(0),medians=double(length(x)),lMedians=as.integer(n))
	median <- res$result
	attr(median, "medians") <- unique(res$medians[1:res$lMedians])
	return(median)
}

medianCircular <- function(x, na.rm=FALSE, type="Fisher", deviation=FALSE, control.circular=list(), ...) {
  .Deprecated(new="median.circular")
  
  ## For now only the definition in
  ## equations 2.32 & 2.33
  ## from N.I. Fisher's 'Statistical Analysis of Circular Data',
  ## Cambridge Univ. Press 1993.
  ## is implemented
   type <- match.arg(type)  
   if (na.rm)
       x <- x[!is.na(x)]
   if (length(x)==0) {
        warning("No observations (at least after removing missing values)")
        return(NULL)
   }

   if (is.circular(x)) {
      datacircularp <- circularp(x)
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
   x <- conversion.circular(x, units="radians")
   attr(x, "class") <- attr(x, "circularp") <-  NULL
  circmedian <- list()
   if (type=="Fisher")
     circmedian$median <- MedianCircularRad(x)
   else
     stop("Others 'type' not yet implemented")
   circmedian$median <- conversion.circular(circular(circmedian$median, template=datacircularp$template, zero=datacircularp$zero, rotation=datacircularp$rotation), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
   if (deviation) {
     circmedian$deviation <- MeanDeviationRad(x)
     return(circmedian)
   }
   else
     return(circmedian$median)
}
