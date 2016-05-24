
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   lm.circular.cl function                                 #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: July, 25, 2006                                    #
#   Version: 0.2-2                                          #
#                                                           #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#############################################################

lm.circular.cl <- function(y, x, init=NULL, verbose=FALSE, tol=1e-10, control.circular=list()) {
    # Handling missing values
    ok <- complete.cases(x, y)
    if (NCOL(x)==1) {
        x <- x[ok]
    } else {
        x <- x[ok,]
    }
    y <- y[ok]
    
    if ((n <- length(y))==0) {
        warning("No observations (at least after removing missing values)")
        return(NULL)
    }
    
    if (is.null(init))
       stop("'init' is missing with no default")
    if (is.vector(x))
       x <- cbind(x)
    if (NCOL(x)!=length(init))
       stop("'init' must have the same number of elements as the columns of 'x'")
    
    if (is.circular(y)) {
       datacircularp <- circularp(y)     
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
    
    y <- conversion.circular(y, units="radians", zero=0, rotation="counter", modulo="2pi")
    attr(y, "circularp") <- attr(y, "class") <- NULL

    result <- LmCircularclRad(y, x, init, verbose, tol)
    result$mu <- conversion.circular(circular(result$mu), dc$units, dc$type, dc$template, dc$modulo, dc$zero, dc$rotation)
    if (dc$units=="degrees") result$se.mu <- result$se.mu*180/pi

    result$call <- match.call()
    class(result) <- "lm.circular.cl"
    return(result)
}

LmCircularclRad <- function(y, x, init, verbose, tol) {
   n <- length(y)
   y <- y%%(2*pi)
   betaPrev <- init  
   S <- sum(sin(y-2*atan(x%*%betaPrev)))/n
   C <- sum(cos(y-2*atan(x%*%betaPrev)))/n
   R <- sqrt(S^2 + C^2)
   mu <- atan2(S,C)
   k  <- A1inv(R)
   diff <- tol+1
   iter <- 0
   while (diff > tol){
      iter <- iter + 1
      u <- k*sin(y - mu - 2*atan(x%*%betaPrev))
      A <- diag(k*A1(k), nrow=n)
      g.p <- diag(apply(x, 1, function(row, betaPrev) 2/(1+(t(betaPrev)%*%row)^2), betaPrev=betaPrev), nrow=n)
      D <- g.p%*%x
      betaNew <- lm(t(D)%*%(u+A%*%D%*%betaPrev) ~ t(D)%*%A%*%D - 1)$coefficients
      diff <- max(abs(betaNew - betaPrev))
      betaPrev <- betaNew
        
      S <- sum(sin(y-2*atan(x%*%betaPrev)))/n
      C <- sum(cos(y-2*atan(x%*%betaPrev)))/n
      R <- sqrt(S^2 + C^2)
      mu <- atan2(S,C)
      k  <- A1inv(R)
        
      if (verbose){
         log.lik <- -n*log(besselI(x = k, nu = 0, expon.scaled = FALSE)) + k*sum(cos(y-mu-2*atan(x%*%betaNew)))
         cat("Iteration ", iter, ":    Log-Likelihood = ", log.lik, "\n")
      }
   }
    
   log.lik <- -n*log(besselI(x = k, nu = 0, expon.scaled = FALSE)) + k*sum(cos(y-mu-2*atan(x%*%betaNew)))
   cov.beta <- solve(t(D)%*%A%*%D)
   se.beta <- sqrt(diag(cov.beta))
   se.kappa <- sqrt(1/(n*(1-A1(k)^2-A1(k)/k)))
   se.mu <- 1/sqrt((n-ncol(x))*k*A1(k))
   t.values <- abs(betaNew/se.beta)
   p.values <- 1-pnorm(t.values)
   betaNew <- as.vector(betaNew)
   result <- list()
   result$x <- x
   result$y <- y
   result$mu <- mu
   result$se.mu <- se.mu
   result$kappa <- k
   result$se.kappa <- se.kappa
   result$coefficients <- betaNew
   result$cov.coef <- cov.beta
   result$se.coef <- se.beta
   result$log.lik <- log.lik
   result$t.values <- t.values
   result$p.values <- p.values
   return(result)
}

#############################################################
#                                                           #
#   print.lm.circular.cl function                           #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: April, 27, 2005                                   #
#   Version: 0.2                                            #
#                                                           #
#   Copyright (C) 2005 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.lm.circular.cl <- function(x, digits = max(3, getOption("digits") - 3), signif.stars= getOption("show.signif.stars"), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="") 
 
    result.matrix <- cbind(x$coefficients, x$se.coef, x$t.values, x$p.values)
    dimnames(result.matrix) <- list(names(x$x),c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    cat("\n Circular-Linear Regression \n")
    cat("\n Coefficients:\n")
    printCoefmat(result.matrix, digits=digits, signif.stars=signif.stars, ...)
    cat("\n")
    cat(" Log-Likelihood: ", format(x$log.lik, digits=digits), "\n")
    cat("\n Summary: (mu in radians)\n")
    cat("  mu: ", format(x$mu, digits=digits), "(", format(x$se.mu, digits=digits), ")  kappa: ",  format(x$kappa, digits=digits), "(", format(x$se.kappa, digits=digits), ")\n")
    cat("p-values are approximated using normal distribution\n\n")
    invisible(x)
}
