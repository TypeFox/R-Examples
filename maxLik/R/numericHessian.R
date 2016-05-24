numericHessian <- function(f, grad=NULL, t0, eps=1e-6, fixed,
                           ...) {
   a <- f(t0, ...)
   if(is.null(grad)) {
      numericNHessian( f = f, t0 = t0, eps = eps, fixed=fixed, ...)
                                        # gradient not provided -> everything numerically
   } else {
      numericGradient( f = grad, t0 = t0, eps = eps, fixed=fixed, ...)
                                        # gradient is provided -> Hessian is grad grad
   }
}

numericNHessian <- function( f, t0, eps=1e-6, fixed, ...) {
   ## Numeric Hessian without gradient
   ## Assume f() returns a scalar
   ## 
   ## fixed   calculate the Hessian only for the non-fixed parameters
   warnMessage <- function(theta, value) {
      ## issue a warning if the function value at theta is not a scalar
      max.print <- 10
      if(length(value) != 1) {
         warnMsg <- "Function value at\n"
         warnMsg <- c(warnMsg,
                      paste(format(theta[seq(length=min(max.print,length(theta)))]),
                             collapse=" "),
                      "\n")
         if(max.print < length(theta))
             warnMsg <- c(warnMsg, "...\n")
         warnMsg <- c(warnMsg, " =\n")
         warnMsg <- c(warnMsg,
                      paste(format(value[seq(length=min(max.print,length(value)))]),
                                   collapse=" "),
                      "\n")
         if(max.print < length(value))
             warnMsg <- c(warnMsg, "...\n")
         warnMsg <- c(warnMsg, "but numeric Hessian only works on numeric scalars\n",
                      "Component set to NA")
         return(warnMsg)
      }
      if(!is.numeric(value))
          stop("The function value must be numeric")
      return(NULL)
   }
   f00 <- f( t0, ...)
   if(!is.null(msg <- warnMessage(t0, f00))) {
      warning(msg)
      f00 <- NA
   }
   eps2 <- eps*eps
   N <- length( t0)
   H <- matrix(NA, N, N)
   if(missing(fixed))
       fixed <- rep(FALSE, length(t0))
   for( i in 1:N) {
      if(fixed[i])
          next
      for( j in 1:N) {
         if(fixed[j])
             next
         t01 <- t0
         t10 <- t0
         t11 <- t0
                                          # initial point
         t01[i] <- t01[i] + eps
         t10[j] <- t10[j] + eps
         t11[i] <- t11[i] + eps
         t11[j] <- t11[j] + eps
         f01 <- f( t01, ...)
         if(!is.null(msg <- warnMessage(t01, f01))) {
            warning(msg)
            f01 <- NA
         }
         f10 <- f( t10, ...)
         if(!is.null(msg <- warnMessage(t10, f10))) {
            warning(msg)
            f10 <- NA
         }
         f11 <- f( t11, ...)
         if(!is.null(msg <- warnMessage(t11, f11))) {
            warning(msg)
            f11 <- NA
         }
         H[i,j] <- ( f11 - f01 - f10 + f00)/eps2
      }
   }
   return( H )
}
