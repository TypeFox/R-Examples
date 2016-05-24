#############################################################
#                                                           #
#   mle.cardioid function                                   #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: December, 6, 2005                                   #
#   Copyright (C) 2005 Claudio Agostinelli                  #
#                                                           #
#   Version 0.1-2                                           #
#############################################################

mle.cardioid <- function(x, mu, rho=0, max.iter=100, tol=1e-3) {

 # Handling missing values
 x <- na.omit(x)
 if (length(x)==0) {
     warning("No observations (at least after removing missing values)")
     return(NULL)
 }
 
 if (rho==0) {
    mu <- NA
    convergence <- FALSE
 } else {
    x <- as.circular(x)
    xcircularp <- circularp(x)
    units <- xcircularp$units
    x <- conversion.circular(x, units="radians")

    n <- length(x)
    sinr <- sum(sin(x))
    cosr <- sum(cos(x))
    if (missing(mu))
        mu <- atan2(sinr, cosr)
    diff <-  tol + 1
    i <- 0
    while (diff>tol & i <= max.iter) {
           i <- i + 1
           mu.old <- mu
           temp <- 1+2*rho*cos(x-mu)
           mu <- atan2(sum(sin(x)/temp),sum(cos(x)/temp))
           cat("i ", i, "\n")
           cat("mu ", mu, "\n")
           cat("temp ", cos(x-mu)[1], "\n")
           diff <- abs(mu-mu.old)
    }

    convergence <- TRUE
    if (i > max.iter)
        convergence <- FALSE
    
    if (units=="degrees") {
        mu <- mu/pi*180
    }
 }
    attr(mu, "circularp") <- xcircularp
    attr(mu, "class") <- "circular"
    result <- list()
    result$call <- match.call()
    result$mu <- mu
    result$convergence <- convergence
    class(result) <- "mle.cardioid"
    return(result)
}

#############################################################
#                                                           #
#	print.mle.cardioid function                         #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: April, 29, 2003                               #
#	Version: 0.1                                        #
#                                                           #
#	Copyright (C) 2003 Claudio Agostinelli              #
#                                                           #
#############################################################

print.mle.cardioid <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("mu: \n")
    print(x$mu, digits=digits)
    cat("\n")
    if (!x$convergence)
        cat("Maximum number of iteration is reached \n")
    invisible(x)
}

ver <- function(x, mu, rho) {
      prod(dcardioid(x, mu, rho))
}

score <- function(x, mu, rho) {
           temp <- 1+2*rho*cos(x-mu)
           cos(mu)*sum(sin(x)/temp)-sin(mu)*sum(cos(x)/temp)
}

#grid <- seq(0, 2*pi, 0.1)
#res <- res.s <- vector(length=0)

#for(i in 1:length(grid)) {
# res <- c(res, ver(x, grid[i], rho))
# res.s <- c(res.s, score(x, grid[i], rho))
#}
