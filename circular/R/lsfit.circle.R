
#############################################################
#                                                           #
#       Original Splus: Ulric Lund                          #
#       E-mail: ulund@calpoly.edu                           #
#                                                           #
#############################################################

#############################################################
#                                                           #
#   lsfit.circle function                                   #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: May, 27, 2006                                     #
#   Version: 0.3                                            #
#                                                           #
#   Copyright (C) 2006 Claudio Agostinelli                  #
#                                                           #
#############################################################

lsfit.circle <- function(x, y=NULL, init=NULL, units=c("radians", "degrees"), template=c("none", "geographics"), modulo=c("asis", "2pi", "pi"), zero=0, rotation=c("counter", "clock"), ...) {
    type <- "angles"
    units <- match.arg(units)
    template <- match.arg(template) 
    modulo <- match.arg(modulo)
    rotation <- match.arg(rotation)

    if (is.null(y)) {
        if (NCOL(x)==2) {
            y <- x[,2]
            x <- x[,1]
        } else {
            stop("Must be either 'x' a matrix with two columns or 'x' and 'y' a vector")
        }
    }
    x <- as.vector(x)
    y <- as.vector(y)
    
    if (length(x)!=length(y)) {
        stop("'x' and 'y' must have the same length") 
    }
    
    # Handling missing values
    ok <- complete.cases(x, y)
    x <- x[ok]
    y <- y[ok]
    if (length(y)==0) {
        warning("No observations (at least after removing missing values)")
        return(NULL)
    }
    
    result <- LsfitCircleRad(x, y, init, ...)
    
    result$angles <- conversion.circular(circular(result$angles), units, type, template, modulo, zero, rotation)
    result$call <- match.call()
    class(result) <- "lsfit.circle"
    return(result)
}


LsfitCircleRad <- function(x, y, init, ...) {
   if (is.null(init)) {
      x.mean <- mean.default(x)
      y.mean <- mean.default(y)
      est.r <- max(c(abs(x-x.mean), abs(y-y.mean)))
      init <- c(est.r, x.mean, y.mean) 
   }
    
   obj.fun <- function(p, x, y){
      sum((p[1]-sqrt((x-p[2])^2 + (y-p[3])^2))^2)
   }
    
   grad.fun <- function(p, x, y){
      n <- length(x); r <- p[1]; a <- p[2]; b <- p[3]
      common <- sqrt((x-a)^2 + (y-b)^2)
      g.e1 <- 2*n*r - 2 * sum(common)
      g.e2 <- 2*n*a - 2*sum(x) + 2 * r * sum((x-a)/common)
      g.e3 <- 2*n*b - 2*sum(y) + 2 * r * sum((y-b)/common)
# NO HESSIAN is used by optim
#        h.e11 <- 2*n
#        h.e21 <- 2*sum((x-a)/common)
#        h.e22 <- 2*n - 2*r*sum((y-b)^2/common^3)
#        h.e31 <- 2*sum((y-b)/common)
#        h.e32 <- 2*r*sum((x-a)*(y-b)/common^3)
#        h.e33 <- 2*n - 2*r*sum((x-a)^2/common^3)
#        pppp <- list(gradient=c(g.e1, g.e2, g.e3), hessian=c(h.e11, h.e21, h.e22, h.e31, h.e32, h.e33))
      return(c(g.e1, g.e2, g.e3))
   } 
#    nlminb(start = init, obj = obj.fun, gradient = grad.fun, hessian=TRUE, x = x, y = y, ...)
   res <- optim(par=init, fn=obj.fun, gr=grad.fun, hessian=TRUE, x = x, y = y, method ="BFGS", ...)
   result <- list()

   result$coefficients <- res$par
   names(result$coefficients) <- c("r", "a", "b") 
   result$x <- x
   result$y <- y
   result$x.centered <- x - res$par[2] 
   result$y.centered <- y - res$par[3]
   result$angles <- atan2(y=result$y.centered, x=result$x.centered)
   result$radius <- sqrt(result$x.centered^2+result$y.centered^2)
   result$convergence <- res$convergence
   result$optim <- res
   return(result) 
}

#############################################################
#                                                           #
#   print.lsfit.circle function                             #
#   Author: Claudio Agostinelli                             #
#   E-mail: claudio@unive.it                                #
#   Date: April, 27, 2005                                   #
#   Version: 0.1-1                                          #
#                                                           #
#   Copyright (C) 2005 Claudio Agostinelli                  #
#                                                           #
#############################################################

print.lsfit.circle <- function(x, digits = max(3, getOption("digits") - 3), ...) {
    cat("\nCall:\n",deparse(x$call),"\n\n",sep="") 
 
    coef <- matrix(x$coefficients, nrow=1)
    dimnames(coef) <- list(" Coeff: ", c("r", "a", "b"))
      print.default(format(coef, digits=digits),
		  print.gap = 2, quote = FALSE)
    cat("\n Summary in Rectangular Coordinates \n")
    print(summary(data.frame(x=x$x, y=x$y)))
    cat("\n Summary in Polar Coordinates of Recentered Observations \n")
    print(summary(data.frame(angles=x$angles, radius=x$radius)))

    if (x$convergence) {
        cat("Warnings: convergence problems in 'optim': ", x$convergence, "\n\n")
    }
    invisible(x)
}

