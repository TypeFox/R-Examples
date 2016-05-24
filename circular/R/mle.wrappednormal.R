#############################################################
#                                                           #
#   mle.wrappednormal function                              #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: August, 10, 2006                                  #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-2                                           #
#############################################################

mle.wrappednormal <- function(x, mu=NULL, rho=NULL, sd=NULL, K=NULL, tol=1e-5, min.sd=1e-3, min.k=10, max.iter=100, verbose=FALSE, control.circular=list()) {

    # Handling missing values
    x <- na.omit(x)
    if (length(x)==0) {
        warning("No observations (at least after removing missing values)")
        return(NULL)
    }
    
    if (is.circular(x)) {
       datacircularp <- circularp(x)     
    } else if (is.circular(mu)) {
              datacircularp <- circularp(mu)     
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
    
    x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
    attr(x, "class") <- attr(x, "circularp") <- NULL
    if (!is.null(mu)) {
       mu <- conversion.circular(mu, units="radians", zero=0, rotation="counter", modulo="2pi")
       attr(mu, "class") <- attr(mu, "circularp") <- NULL
    }

    res <- MlewrappednormalRad(x, mu, rho, sd, min.sd, K, min.k, tol, max.iter, verbose)

    mu <- conversion.circular(circular(res[1]), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)

    result <- list()
    result$call <- match.call()
    result$mu <- mu
    result$rho <- res[2]
    result$sd <- res[3]
    result$est.mu <- res[4]
    result$est.rho <- res[5]
    result$convergence <- TRUE
    if (res[6] > max.iter) {
        result$convergence <- FALSE
    }
    class(result) <- "mle.wrappednormal"
    return(result)
}

MlewrappednormalRad <- function(x, mu=NULL, rho=NULL, sd=NULL, min.sd, K=NULL, min.k=10, tol, max.iter, verbose) {
   n <- length(x)
   sinr <- sum(sin(x))
   cosr <- sum(cos(x))

   est.mu <- FALSE 
   if (is.null(mu)) {  
       mu <- atan2(sinr, cosr)
       est.mu <- TRUE
   }
   est.rho <- FALSE
   if (is.null(sd)) {
      if (is.null(rho)) {
         sd <- sqrt(-2*log(sqrt(sinr^2 + cosr^2)/n))
         if (is.na(sd) || sd < min.sd)
            sd <- min.sd
         est.rho <- TRUE
      } else {
         sd <- sqrt(-2*log(rho))
      }
   }
     
   xdiff <- 1+tol
   iter <- 0
   if (is.null(K)) {
      range <- max(mu, x) - min(mu, x)
      K <- (range+6*sd)%/%(2*pi)+1
      K <- max(min.k, K)
   }
    
   while (xdiff > tol & iter <= max.iter) {
      iter <- iter + 1
      mu.old <- mu
      sd.old <- sd
      z <- .Fortran("mlewrpno",
         as.double(x),
         as.double(mu),
         as.double(sd),
         as.integer(n),
         as.integer(K),
         as.integer(est.mu),
         as.integer(est.rho),
         w=double(n),
         wk=double(n),
         wm=double(n),
         PACKAGE="circular"
      )
      w <- z$w
      wk <- z$wk
      wm <- z$wm
           
      if (est.mu) {
         mu <- sum(x)/n
         if (any(wk!=0)) {
            mu <- mu + 2*pi*mean.default(wk[wk!=0]/w[wk!=0])
         }
      }
      if (est.rho) {
         if (any(wm!=0)) {
            sd <- sqrt(sum(wm[wm!=0]/w[wm!=0])/n)
         } else {
            sd <- min.sd
         }
      }

      if (verbose) {
         cat("mu: ", mu, "\n")
         cat("rho: ", exp(-sd^2/2), "\n")              
         cat("sd: ", sd, "\n")
      }
      xdiff <- max(abs(mu - mu.old), abs(sd - sd.old))
  }

  rho <- exp(-sd^2/2)
  result <- c(mu, rho, sd, est.mu, est.rho, iter)
  return(result)
}

#############################################################
#                                                           #
#	print.mle.wrappednormal function                    #
#	Author: Claudio Agostinelli                         #
#	E-mail: claudio@unive.it                            #
#	Date: November, 19, 2003                                #
#	Version: 0.1-2                                      #
#                                                           #
#	Copyright (C) 2003 Claudio Agostinelli              #
#                                                           #
#############################################################

print.mle.wrappednormal <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("mu: ")
    cat(format(x$mu, digits=digits), "\n")
    cat("\n")
    cat("rho: ")    
    cat(format(x$rho, digits=digits), "\n")
    cat("\n")
    cat("sd: ")       
    cat(format(x$sd, digits=digits), "\n")
    cat("\n")   
    if (!x$est.mu) cat("mu is known\n")
    if (!x$est.rho) {
        cat("rho and sd are known\n")
    }
    if (!x$convergence) cat("\nThe convergence is not achieved after the prescribed number of iterations \n")
    invisible(x)
}

