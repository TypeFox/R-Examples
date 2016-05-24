#############################################################
#                                                           #
#	wle.weights function                                    #
#	Author: Claudio Agostinelli                             #
#	E-mail: claudio@unive.it                                #
#	Date: April, 02, 2002                                   #
#	Version: 0.1-4                                          #
#                                                           #
#	Copyright (C) 2002 Claudio Agostinelli                  #
#                                                           #
#############################################################

wle.weights <- function(x, y=NULL, smooth=0.0031, sigma2, raf=1, location=FALSE, max.iter=1000, tol=10^(-6)) {

    result <- list()

    if (is.null(y)) {
        y <- x
    } else {
        location <- FALSE
    }  

    loc <- 0
    loc.old <- loc + tol + 1
    nx <- length(x)
    ny <- length(y)
    iter <- 0
    conv <- TRUE

    while (abs(loc-loc.old)>tol & conv) {

           iter <- iter + 1
           loc.old <- loc
           if (location) {
               y <- z <- x - loc
           } else {
               z <- x
           }

           w.temp <- .Fortran("wlew",
	              as.double(z), 
	              as.integer(nx),
	              as.double(y), 
	              as.integer(ny), 
	              as.integer(raf),
	              as.double(smooth),
	              as.double(sigma2),
	              totweight=double(1),
	              weights=double(nx),
				  PACKAGE="wle")

           loc <- sum(w.temp$weights*x)/sum(w.temp$weights)

           if (!location) loc <- loc.old
           if (iter > max.iter) conv <- FALSE
    }

    result$weights <- w.temp$weights
    if (location) {
        result$location <- loc
    } else {
        result$location <- NA
    }
    result$conv <- conv

return(result)
}

