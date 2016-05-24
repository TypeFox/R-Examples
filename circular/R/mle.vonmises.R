
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   mle.vonmises function                                   #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: August, 10, 2006                                  #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.3-2                                           #
#############################################################

mle.vonmises <- function(x, mu=NULL, kappa=NULL, bias=FALSE, control.circular=list()) {

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
    
    res <- MlevonmisesRad(x, mu, kappa, bias)

    mu <- conversion.circular(circular(res[1]), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
    if (dc$units=="degrees") res[2] <- res[2]*180/pi

    result <- list()
    result$call <- match.call()
    result$mu <- mu
    result$kappa <- res[4]
    result$se.mu <- res[2]
    result$se.kappa <- res[5]
    result$est.mu <- res[3]
    result$est.kappa <- res[6]
    result$bias <- bias
    class(result) <- "mle.vonmises"
    return(result)
}

MlevonmisesRad <- function(x, mu=NULL, kappa=NULL, bias=FALSE) {
    n <- length(x)
    sinr <- sum(sin(x))
    cosr <- sum(cos(x))
    est.mu <- FALSE 
    if (is.null(mu)) {  
        mu <- atan2(sinr, cosr)
        est.mu <- TRUE
    }
    est.kappa <- FALSE
    if (is.null(kappa)) {
        V <- mean.default(cos(x - mu))
        if (V > 0) {
            kappa <- A1inv(V)
        } else {
            kappa <- 0
        }
        if (bias == TRUE) {
            if (kappa < 2) {
                kappa <- max(kappa - 2 * (n * kappa)^-1, 0)
            } else {
                kappa <- ((n - 1)^3 * kappa)/(n^3 + n)
            }
        }
        est.kappa <- TRUE
    }

    A1temp <- A1(kappa)
    se.mu <- se.kappa <- 0
    if (est.mu) se.mu <- sqrt(1/(n*kappa*A1temp))
    if (est.kappa) se.kappa <- sqrt(1/(n*(1-A1temp/kappa-A1temp^2)))
    result <- c(mu, se.mu, est.mu, kappa, se.kappa, est.kappa)
    return(result)
}

#############################################################
#                                                           #
#   print.mle.vonmises function                             #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: May, 22, 2006                                     #
#   Version: 0.2                                            #
#                                                           #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.mle.vonmises <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("mu: ")
    cat(format(x$mu, digits=digits), " (", format(x$se.mu, digits=digits), ")\n")
    cat("\n")
    cat("kappa: ")    
    cat(format(x$kappa, digits=digits), " (", format(x$se.kappa, digits=digits), ")\n")
    cat("\n")    
    if (!x$est.mu) cat("mu is known\n")
    if (!x$est.kappa) cat("kappa is known\n")
    if (x$bias) cat("Bias correction (Best and Fisher, 1981) applied to kappa\n")

    invisible(x)
}
