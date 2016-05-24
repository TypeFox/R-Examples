#####################################################################
#                                                                   #
#   trigonometric.polynomials function                              #
#   Author: Claudio Agostinelli and Alessandro Gagliardi            #
#   Email: claudio@unive.it                                         #
#   Date: January, 04, 2013                                         #
#   Copyright (C) 2013 Claudio Agostinelli and Alessandro Gagliardi #
#                                                                   #
#   Version 0.1                                                     #
#####################################################################

trigonometric.polynomials <- function(x, p = 1, center = FALSE) {
  p <- as.vector(p)
  x <- conversion.circular(x, units="radians", zero=0, rotation="counter")
  attr(x, "class") <- attr(x, "circularp") <-  NULL
  x <- as.matrix(x)
  if (is.null(colnames(x)))
    colnam <- paste('x', 1L:NCOL(x), sep='')
  else
    colnam <- colnames(x)
#  if (is.null(colnames(x))) {
#    if (NCOL(x)==1)
#      colnames(x) <- deparse(substitute(x))
#    else
#      colnames(x) <- paste(deparse(substitute(x)), 1L:NCOL(x), sep='')
#  }
  result <- matrix(NA, nrow=NROW(x), ncol=0)
  for (i in 1L:NCOL(x)) {
    for (j in 1L:length(p)) {
      res <- TrigonometricPolynomialsRad(x[,i], p[j], center)
      colnames(res) <- c(paste('cos(', ifelse(p[j]==1,'',round(p[j],3)), colnam[i], ifelse(center,'-mean',''), ')', sep=''), paste('sin(', ifelse(p[j]==1,'',round(p[j],3)), colnam[i], ifelse(center,'-mean',''), ')', sep=''))
      result <- cbind(result, res)
    }
  }
  return(result)
}

TP <- function(x, p = 1, center = FALSE) {
  tp <- trigonometric.polynomials(x = x, p = p, center = center)
  class(tp) <- unique(c("AsIs", oldClass(tp)))
  return(tp)
}
  
TrigonometricPolynomialsRad <- function(x, p, center) {
  center <- as.numeric(center)
  sinr <- sum(sin(x))
  cosr <- sum(cos(x))
  circmean <- atan2(sinr, cosr)
  sin.p <- sin(p * (x - circmean * center))
  cos.p <- cos(p * (x - circmean * center))
  result <- cbind(cos.p, sin.p)
  return(result)
}

