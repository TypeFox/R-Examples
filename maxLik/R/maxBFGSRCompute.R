maxBFGSRCompute <- function(fn,
                         start, 
                            finalHessian=TRUE,
                  fixed=NULL,
                            control=maxControl(),
                  ...) {
   ## This function is originally developed by Yves Croissant (and placed in 'mlogit' package).
   ## Fitted for 'maxLik' by Ott Toomet, and revised by Arne Henningsen
   ## 
   ## BFGS maximisation, implemented by Yves Croissant
   ## Parameters:
   ## fn          - the function to be minimized.  Returns either scalar or
   ##               vector value with possible attributes 
   ##               constPar and newVal
   ##               fn must return the value with attribute 'gradient'
   ##               (and also attribute 'hessian' if it should be returned)
   ##               fn must have an argument sumObs
   ## start       - initial parameter vector (eventually w/names)
   ## finalHessian  include final Hessian?  As computing final hessian does not carry any extra penalty for NR method, this option is
   ##               mostly for compatibility reasons with other maxXXX functions.
   ##               TRUE/something else  include
   ##               FALSE                do not include
   ## fixed       - a logical vector -- which parameters are taken as fixed.
   ## control       MaxControl object:
   ##     steptol     - minimum step size
   ##     lambdatol   - max lowest eigenvalue when forcing pos. definite H
   ##     qrtol       - tolerance for qr decomposition
   ##     qac           How to handle the case where new function value is
   ##               smaller than the original one:
   ##                  "stephalving"   smaller step in the same direction
   ##                  "marquardt"     Marquardt (1963) approach
   ##     The stopping criteria
   ##     tol         - maximum allowed absolute difference between sequential values
   ##     reltol      - maximum allowed reltive difference (stops if < reltol*(abs(fn) + reltol)
   ##     gradtol     - maximum allowed norm of gradient vector
   ## 
   ##     iterlim     - maximum # of iterations
   ##     
   ##               Other paramters are treated as variable (free).
   ##
   ## RESULTS:
   ## a list of class "maxim":
   ## maximum     function value at maximum
   ## estimate    the parameter value at maximum
   ## gradient        gradient
   ## hessian         Hessian
   ## code        integer code of success:
   ##             1 - gradient close to zero
   ##             2 - successive values within tolerance limit
   ##             3 - could not find a higher point (step error)
   ##             4 - iteration limit exceeded
   ##             100 - initial value out of range
   ## message     character message describing the code
   ## last.step   only present if code == 3 (step error).  A list with following components:
   ##             theta0    - parameter value which led to the error
   ##             f0        - function value at these parameter values
   ##             climb     - the difference between theta0 and the new approximated parameter value (theta1)
   ##             fixed     - logical vector, which parameters are constant (fixed, inactive, non-free)
   ## fixed       logical vector, which parameters were treated as constant (fixed, inactive, non-free)
   ## iterations  number of iterations
   ## type        "BFGSR maximisation"
   ## 

   ##
   max.eigen <- function( M) {
      ## return maximal eigenvalue of (symmetric) matrix
      val <- eigen(M, symmetric=TRUE, only.values=TRUE)$values
      val[1]
      ## L - eigenvalues in decreasing order, [1] - biggest in abs value
   }
   ##
   maxim.type <- "BFGSR maximization"
  param <- start
   nimed <- names(start)
   nParam <- length(param)
   ##
  chi2 <- 1E+10
  iter <- 0
  # eval a first time the function, the gradient and the hessian
  x <- sumKeepAttr( fn( param, fixed = fixed, sumObs = FALSE,
      returnHessian = FALSE, ... ) )
            # sum of log-likelihood value but not sum of gradients
   if (slot(control, "printLevel") > 0)
    cat( "Initial value of the function :", x, "\n" )
   if(is.na(x)) {
      result <- list(code=100, message=maximMessage("100"),
                     iterations=0,
                     type=maxim.type)
      class(result) <- "maxim"
      return(result)
   }
   if(is.infinite(x) & (x > 0)) {
                                        # we stop at +Inf but not at -Inf
      result <- list(code=5, message=maximMessage("5"),
                     iterations=0,
                     type=maxim.type)
      class(result) <- "maxim"
      return(result)
   }
   if( isTRUE( attr( x, "gradBoth" ) ) ) {
      warning( "the gradient is provided both as attribute 'gradient' and",
         " as argument 'grad': ignoring argument 'grad'" )
   }
   if( isTRUE( attr( x, "hessBoth" ) ) ) {
      warning( "the Hessian is provided both as attribute 'hessian' and",
         " as argument 'hess': ignoring argument 'hess'" )
   }
   ##
   ## gradient by individual observations, used for BHHH approximation of initial Hessian.
   ## If not supplied by observations, we use the summed gradient.
   gri <- attr( x, "gradient" )
   gr <- sumGradients( gri, nParam = length( param ) )
   if(slot(control, "printLevel") > 2) {
      cat("Initial gradient value:\n")
      print(gr)
   }
   if(any(is.na(gr[!fixed]))) {
      stop("NA in the initial gradient")
   }
   if(any(is.infinite(gr[!fixed]))) {
      stop("Infinite initial gradient")
   }
   if(length(gr) != nParam) {
      stop( "length of gradient (", length(gr),
         ") not equal to the no. of parameters (", nParam, ")" )
   }
   ## initial approximation for inverse Hessian.  We only work with the non-fixed part
   if(observationGradient(gri, length(param))) {
      invHess <- -solve(crossprod(gri[,!fixed]))
                           # initial approximation of inverse Hessian (as in BHHH), if possible
      if(slot(control, "printLevel") > 3) {
         cat("Initial inverse Hessian by gradient crossproduct\n")
         if(slot(control, "printLevel") > 4) {
            print(invHess)
         }
      }
   }
   else {
      invHess <- -1e-5*diag(1, nrow=length(gr[!fixed]))
                           # ... if not possible (Is this OK?).  Note we make this negative definite.
      if(slot(control, "printLevel") > 3) {
         cat("Initial inverse Hessian is diagonal\n")
         if(slot(control, "printLevel") > 4) {
            print(invHess)
         }
      }
   }
   if( slot(control, "printLevel") > 1) {
      cat("-------- Initial parameters: -------\n")
      cat( "fcn value:",
      as.vector(x), "\n")
      a <- cbind(start, gr, as.integer(!fixed))
      dimnames(a) <- list(nimed, c("parameter", "initial gradient",
                                          "free"))
      print(a)
      cat("------------------------------------\n")
   }
   samm <- NULL
                           # this will be returned in case of step getting too small
   I <- diag(nParam - sum(fixed))
  direction <- rep(0, nParam)
   ## ----------- Main loop ---------------
   repeat {
      iter <- iter + 1
      if( iter > slot(control, "iterlim")) {
         code <- 4; break
      }
      if(any(is.na(invHess))) {
         cat("Error in the approximated (free) inverse Hessian:\n")
         print(invHess)
         stop("NA in Hessian")
      }
      if(slot(control, "printLevel") > 0) {
         cat("Iteration ", iter, "\n")
         if(slot(control, "printLevel") > 3) {
            cat("Eigenvalues of approximated inverse Hessian:\n")         
            print(eigen(invHess, only.values=TRUE)$values)
            if(slot(control, "printLevel") > 4) {
               cat("inverse Hessian:\n")
               print(invHess)
            }
         }
      }
      ## Next, ensure that the approximated inverse Hessian is negative definite for computing
      ## the new climbing direction.  However, retain the original, potentially not negative definite
      ## for computing the following approximation.
      ## This procedure seems to work, but unfortunately I have little idea what I am doing :-(
      approxHess <- invHess
                           # approxHess is used for computing climbing direction, invHess for next approximation
      while((me <- max.eigen( approxHess)) >= -slot(control, "lambdatol") |
         (qRank <- qr(approxHess, tol=slot(control, "qrtol"))$rank) < sum(!fixed)) {
                                        # maximum eigenvalue -> negative definite
                                        # qr()$rank -> singularity
         lambda <- abs(me) + slot(control, "lambdatol") + min(abs(diag(approxHess)))/1e7
                           # The third term corrects numeric singularity.  If diag(H) only contains
                           # large values, (H - (a small number)*I) == H because of finite precision
         approxHess <- approxHess - lambda*I
         if(slot(control, "printLevel") > 4) {
            cat("Not negative definite.  Subtracting", lambda, "* I\n")
            cat("Eigenvalues of new approximation:\n")         
            print(eigen(approxHess, only.values=TRUE)$values)
            if(slot(control, "printLevel") > 5) {
               cat("new Hessian approximation:\n")
               print(approxHess)
            }
         }
                           # how to make it better?
      }
      ## next, take a step of suitable length to the suggested direction
      step <- 1
      direction[!fixed] <- as.vector(approxHess %*% gr[!fixed])
    oldx <- x
     oldgr <- gr
    oldparam <- param
    param[!fixed] <- oldparam[!fixed] - step * direction[!fixed]
    x <- sumKeepAttr( fn( param, fixed = fixed, sumObs = FALSE,
      returnHessian = FALSE, ... ) )
                           # sum of log-likelihood value but not sum of gradients
      ## did we end up with a larger value?
      while((is.na(x) | x < oldx) & step > slot(control, "steptol")) {
         step <- step/2
         if(slot(control, "printLevel") > 2) {
            cat("Function decreased. Function values: old ", oldx, ", new ", x, ", difference ", x - oldx, "\n")
            if(slot(control, "printLevel") > 3) {
               resdet <- cbind(param = param, gradient = gr, direction=direction, active=!fixed)
               cat("Attempted parameters:\n")
               print(resdet)
            }
            cat(" -> step ", step, "\n", sep="")
         }
       param[!fixed] <- oldparam[!fixed] - step * direction[!fixed]
       x <- sumKeepAttr( fn( param, fixed = fixed, sumObs = FALSE,
         returnHessian = FALSE, ... ) )
            # sum of log-likelihood value but not sum of gradients
    }
    if(step < slot(control, "steptol")) {
                           # we did not find a better place to go...
       samm <- list(theta0=oldparam, f0=oldx, climb=direction)
    }
     gri <- attr( x, "gradient" )
                           # observation-wise gradient.  We only need it in order to compute the BHHH Hessian, if asked so.
     gr <- sumGradients( gri, nParam = length( param ) )
      incr <- step * direction
      y <- gr - oldgr
      if(all(y == 0)) {
                           # gradient did not change -> cannot proceed
         code <- 9; break
      }
      ## Compute new approximation for the inverse hessian
      update <- outer( incr[!fixed], incr[!fixed]) *
          (sum(y[!fixed] * incr[!fixed]) +
           as.vector( t(y[!fixed]) %*% invHess %*% y[!fixed])) / sum(incr[!fixed] * y[!fixed])^2 +
               (invHess %*% outer(y[!fixed], incr[!fixed])
                + outer(incr[!fixed], y[!fixed]) %*% invHess)/
                    sum(incr[!fixed] * y[!fixed])
      invHess <- invHess - update
      ##
      chi2 <- -  crossprod(direction[!fixed], oldgr[!fixed])
    if (slot(control, "printLevel") > 0){
       cat("step = ",step, ", lnL = ", x,", chi2 = ",
           chi2, ", function increment = ", x - oldx, "\n",sep="")
       if (slot(control, "printLevel") > 1){
          resdet <- cbind(param = param, gradient = gr, direction=direction, active=!fixed)
          print(resdet)
          cat("--------------------------------------------\n")
       }
    }
      if( step < slot(control, "steptol")) {
         code <- 3; break
      }
      if( sqrt( crossprod( gr[!fixed] ) ) < slot(control, "gradtol") ) {
         code <-1; break
      }
      if(x - oldx < slot(control, "tol")) {
         code <- 2; break
      }
      if(x - oldx < slot(control, "reltol")*(x + slot(control, "reltol"))) {
         code <- 8; break
      }
      if(is.infinite(x) & x > 0) {
         code <- 5; break
      }
   }
   if( slot(control, "printLevel") > 0) {
      cat( "--------------\n")
      cat( maximMessage( code), "\n")
      cat( iter, " iterations\n")
      cat( "estimate:", param, "\n")
      cat( "Function value:", x, "\n")
   }
   if( is.matrix( gr ) ) {
      if( dim( gr )[ 1 ] == 1 ) {
         gr <- gr[ 1, ]
      }
   }
   names(gr) <- names(param)
   # calculate (final) Hessian
   if(tolower(finalHessian) == "bhhh") {
      if(observationGradient(gri, length(param))) {
          hessian <- - crossprod( gri )
          attr(hessian, "type") <- "BHHH"
      } else {
         hessian <- NULL
         warning("For computing the final Hessian by 'BHHH' method, the log-likelihood or gradient must be supplied by observations")
      }
   } else if(finalHessian) {
      hessian <- attr( fn( param, fixed = fixed, returnHessian = TRUE, ... ) ,
         "hessian" )
   } else {
       hessian <- NULL
   }
   if( !is.null( hessian ) ) {
      rownames( hessian ) <- colnames( hessian ) <- nimed
   }
   ## remove attributes from final value of objective (likelihood) function
   attributes( x )$gradient <- NULL
   attributes( x )$hessian <- NULL
   attributes( x )$gradBoth <- NULL
   attributes( x )$hessBoth <- NULL
   ##
   result <-list(
                  maximum = unname( drop( x ) ),
                  estimate=param,
                  gradient=gr,
                 hessian=hessian,
                  code=code,
                  message=maximMessage( code),
                  last.step=samm,
                                        # only when could not find a
                                        # lower point
                  fixed=fixed,
                  iterations=iter,
                  type=maxim.type)
   if(observationGradient(gri, length(param))) {
       colnames( gri ) <- names( param )
       result$gradientObs <- gri
   }
   result <- c(result, control=control)
                           # attach the control parameters
   class(result) <- c("maxim", class(result))
   invisible(result)
}
