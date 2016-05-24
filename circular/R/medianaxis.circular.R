medianaxis <- function(x, na.rm, ...) UseMethod("medianaxis")
medianaxis.default <- function(x, na.rm, ...) .NotYetImplemented()

#############################################################
#                                                           
#   medianaxis.circular function                                  
#   Author: Claudio Agostinelli and Alessandro Gagliardi
#   E-mail: claudio@unive.it                                
#   Date: August, 3, 2011
#   Version: 0.1                                          
#                                                           
#   Copyright (C) 2012 Claudio Agostinelli and Alessandro Gagliardi
#                                                           
#############################################################

medianaxis.circular <- function(x, na.rm=FALSE, ...) {
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
  circmedianaxis <- MedianAxisCircularRad(x)
  circmedianaxis <- conversion.circular(circular(circmedianaxis$medianaxis), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
  return(circmedianaxis)
}

MedianAxisCircularRad <- function(x) {
  n <- length(x)
  y <- c(x, x+pi)
  y <- sort(y %% (2*pi))
  if (odd <- as.logical(n%%2)) {
  ## odd number of observations
    ismedianaxis <- rep(FALSE, n)
    for (i in 1:length(x)) {
      z <- MinusPiPlusPiRad((x-x[i])%%(2*pi))
      pos <- sum(z > 0 & z != pi)
      neg <- sum(z < 0 & z !=-pi)
      zero <- sum(z==0)
      atpi <- sum(z==pi | z==-pi)
###      ismedianaxis[i] <- pos == (neg)
    }
  } else {
  ## even number of observations    
    eps <- min(diff(y), y[1]+(2*pi-y[2*n]))/4 


  }
}


