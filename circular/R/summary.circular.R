
###############################################################
#                                                             #
#       Original Splus: Ulric Lund                            #
#       E-mail: ulund@calpoly.edu                             #
#                                                             #
###############################################################

#############################################################
#                                                           #
#   summary.circular                                        #
#   Authors: Claudio Agostinelli, David Andel,					#
#   Alessandro Gagliardi                                    #
#   Email: claudio@unive.it, andel@ifi.unizh.ch             #
#   Date: February, 3, 2013                                 #
#   Copyright (C) 2003 Claudio Agostinelli, David Andel     #
#   Copyright (C) 2013 Claudio Agostinelli     					#
#                                                           #
#   Version 0.4                                           	#
#############################################################

summary.circular <- function(object, digits = max(3, getOption("digits") -  3), ...) {
  if (is.matrix(object)) {
    return(summary.matrix(object, ...))
  }
  if (is.data.frame(object)) {
    return(summary.data.frame(object, ...))
  } else {
    nas <- is.na(object)
    object <- object[!nas]
    n <- length(object)
	 qq <- minusPiPlusPi(quantile.circular(object))
	 qq <- signif(c(n, qq[1L:3L], mean.circular(object), qq[4L:5L], rho.circular(object)), digits)
    names(qq) <- c("n", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.",  "Rho")
    if(any(nas))
      c(qq, "NA's" = sum(nas))
    else qq
  }
}

