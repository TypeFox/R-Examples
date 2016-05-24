numericGradient <- function(f, t0, eps=1e-6, fixed,
                            ...) {
   ## numeric gradient of a vector-valued function
   ## f           function, return Nval x 1 vector of values
   ## t0          NPar x 1 vector of parameters
   ## fixed       calculate the gradient based on these parameters only
   ## return:
   ## NvalxNPar matrix, gradient
   ## gradient along parameters which are not active are NA
   warnMessage <- function(theta, value, i) {
      ## issue a warning if the function value at theta is not a scalar
      max.print <- 10
      if(length(value) != nVal) {
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
         warnMsg <- c(warnMsg, "(length ", length(value), ") does not conform with ",
                      "the length at original value ", nVal, "\n")
         warnMsg <- c(warnMsg, "Component ", i, " set to NA")
         return(warnMsg)
      }
      if(!all(is.na(value)) & !is.numeric(value))
          stop("The function value must be numeric for 'numericGradient'")
      return(NULL)
   }
   NPar <- length(t0)
   nVal <- length(f0 <- f(t0, ...))
   grad <- matrix(NA, nVal, NPar)
   row.names(grad) <- names(f0)
   colnames(grad) <- names(t0)
   if(missing(fixed))
       fixed <- rep(FALSE, NPar)
   for(i in 1:NPar) {
      if(fixed[i])
          next
      t2 <- t1 <- t0
      t1[i] <- t0[i] - eps/2
      t2[i] <- t0[i] + eps/2
      ft1 <- f(t1, ...)
      ft2 <- f(t2, ...)
      ## give meaningful error message if the functions give vectors
      ## of different length at t1, t2
      if(!is.null(msg <- warnMessage(t1, ft1, i))) {
         warning(msg)
         ft1 <- NA
      }
      if(!is.null(msg <- warnMessage(t2, ft2, i))) {
         warning(msg)
         ft2 <- NA
      }
      grad[,i] <- (ft2 - ft1)/eps
   }
   return(grad)
}

