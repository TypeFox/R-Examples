
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   mle.wrappedcauchy function                              #
#   Author: Claudio Agostinelli                             #
#   Email: claudio@unive.it                                 #
#   Date: August, 10, 2006                                  #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#   Version 0.2-2                                           #
#############################################################

mle.wrappedcauchy <- function(x, mu=NULL, rho=NULL, tol = 1e-015, max.iter = 100, control.circular=list()) {

    if (tol <= 0) stop("'tol' must be positive")
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

    res <- MlewrappedcauchyRad(x, mu, rho, tol, max.iter)

    mu <- conversion.circular(circular(res[1]), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
    
    result <- list()
    result$call <- match.call()
    result$mu <- mu 
    result$rho <- res[2]
    result$est.mu <- res[3]
    result$est.rho <- res[4]
    result$convergence <- TRUE
    if (!is.na(res[5]) && res[5] > max.iter) {
        result$convergence <- FALSE
    }
    class(result) <- "mle.wrappedcauchy"
    return(result)
}

MlewrappedcauchyRad <- function(x, mu, rho, tol, max.iter) {
    n <- length(x)
    est.mu <- FALSE
    if (is.null(mu)) {
        mu <- MeanCircularRad(x)
        est.mu <- TRUE
    }
    est.rho <- FALSE
    if (is.null(rho)) {
       rho <- RhoCircularRad(x)
       est.rho <- TRUE
    }
    if (rho < 0 | rho > 1) stop("'rho' must be between 0 and 1")
    if (est.mu) {
       if (est.rho) {
          mu1.old <- (2 * rho * cos(mu))/(1 + rho^2)
          mu2.old <- (2 * rho * sin(mu))/(1 + rho^2)
          w.old <- 1/(1 - mu1.old * cos(x) - mu2.old * sin(x))
          flag <- TRUE
          iter <- 0
          while (flag & iter <= max.iter) {
             iter <- iter + 1
             mu1.new <- sum(w.old * cos(x))/sum(w.old)
             mu2.new <- sum(w.old * sin(x))/sum(w.old)
             diff1 <- abs(mu1.new - mu1.old)
             diff2 <- abs(mu2.new - mu2.old)
             if ((diff1 < tol) && (diff2 < tol))
                flag <- FALSE
             else {
                mu1.old <- mu1.new
                mu2.old <- mu2.new
                w.old <- 1/(1 - mu1.old * cos(x) - mu2.old * sin(x))
             }
          }
          mu.const <- sqrt(mu1.new^2 + mu2.new^2)
          mu <- atan2(mu2.new, mu1.new) %% (2 * pi)
          rho <- (1 - sqrt(1 - mu.const^2))/mu.const
       } else {
          score <- function(x, data, rho) {
             sum(sin(data-x)/(1+rho^2-2*rho*cos(data-x)))
          }
          res <- uniroot(f=score, lower=mu-pi/2, upper=mu+pi/2, data=x, rho=rho, tol=tol)
          mu <- res$root
          iter <- NA
       }    
     } else {
         if (est.rho) {
             wt <- function(x, mu, rho) {
                ((1-rho^2)*(1+rho^2-2*rho*cos(x-mu)))^(-1)
             }
             diff <- 1+tol
             iter <- 0
             rho.old <- rho
             while (diff >= tol & iter <= max.iter) {
                iter <- iter + 1
                w <- wt(x, mu, rho)
                sumw <- sum(w)
                sumwcos <- w%*%cos(x-mu)
                rho <- (sumw - sqrt(sumw^2 - sumwcos^2))/sumwcos
                diff <- abs(rho - rho.old)
#####                cat("iter: ", iter, " rho: ", rho, "\n")
                rho.old <- rho
             }
           
         }
     }

    result <- c(mu, rho, est.mu, est.rho, iter)
    return(result)
}

#############################################################
#                                                           #
#   print.mle.wrappednormal function                        #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: May, 22, 2006                                     #
#   Version: 0.2                                            #
#                                                           #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.mle.wrappedcauchy <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
    cat("mu: ")
    cat(format(x$mu, digits=digits), "\n")
    cat("\n")
    cat("rho: ")    
    cat(format(x$rho, digits=digits), "\n")
    cat("\n")    
    if (!x$est.mu) cat("mu is known\n")
    if (!x$est.rho) cat("rho is known\n")
    if (!x$convergence) cat("\n The convergence is not achieved after the prescribed number of iterations \n")
    invisible(x)
}
