maxNRCompute <- function(fn,
                         start, 
                           # maximum lambda for Marquardt (1963)
                         finalHessian=TRUE,
                  bhhhHessian = FALSE,
                  fixed=NULL,
                         control=maxControl(),
                  ...) {
   ## Newton-Raphson maximisation
   ## Parameters:
   ## fn          - the function to be maximized.  Returns either scalar or
   ##               vector value with possible attributes 
   ##               constPar and newVal
   ##               fn must return the value with attributes 'gradient'
   ##               and 'hessian'
   ##               fn must have an argument sumObs
   ## start       - initial parameter vector (eventually w/names)
   ## control       MaxControl object:
   ##     steptol     - minimum step size
   ##     lambda0       initial Hessian corrector (see Marquardt, 1963, p 438)
   ##     lambdaStep    how much Hessian corrector lambda is changed between
   ##                   two lambda trials
   ##                  (nu in Marquardt (1963, p 438)
   ##     maxLambda     largest possible lambda (if exceeded will give step error)
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
   ## finalHessian  include final Hessian?  As computing final hessian does not carry any extra penalty for NR method, this option is
   ##               mostly for compatibility reasons with other maxXXX functions.
   ##               TRUE/something else  include
   ##               FALSE                do not include
   ## fixed       - a logical vector -- which parameters are taken as fixed.
   ##               Other paramters are treated as variable (free).
   ## ...           additional argument to 'fn'.  This may include
   ##               'fnOrig', 'gradOrig', 'hessOrig' if called fromm
   ##               'maxNR'.
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
   ##             5 - infinite function value
   ##             6  infinite gradient
   ##             7  infinite Hessian
   ##             100 - initial value out of range
   ## message     character message describing the code
   ## last.step   only present if code == 3 (step error).  A list with following components:
   ##             theta0    - parameter value which led to the error
   ##             f0        - function value at these parameter values
   ##             climb     - the difference between theta0 and the new approximated parameter value (theta1)
   ##             fixed     - logical vector, which parameters are constant (fixed, inactive, non-free)
   ## fixed       logical vector, which parameters were treated as constant (fixed, inactive, non-free)
   ## iterations  number of iterations
   ## type        "Newton-Raphson maximisation"
   ##
   ## References:
   ## Marquardt (1963), "An algorithm for least-squares estimation of nonlinear
   ##      parameters", J. Soc. Indust. Appl. Math 11(2), 431-441
   ##      
   max.eigen <- function( M) {
      ## return maximal eigenvalue of (symmetric) matrix
      val <- eigen(M, symmetric=TRUE, only.values=TRUE)$values
      val[1]
      ## L - eigenvalues in decreasing order, [1] - biggest in abs value
   }
   ## -------------------------------------------------
   if(slot(control, "qac") == "marquardt")
      marquardt <- TRUE
   else
      marquardt <- FALSE
   ##
   maximType <- "Newton-Raphson maximisation"
   if(marquardt) {
      maximType <- paste(maximType, "with Marquardt (1963) Hessian correction")
   }
   nimed <- names(start)
   nParam <- length(start)
   samm <- NULL
                           # data for the last step that could not find a better
                           # value
   I <- diag(rep(1, nParam))
                           # I is unit matrix
   start1 <- start
   iter <- 0
   returnHessian <- ifelse( bhhhHessian, "BHHH", TRUE )
   f1 <- fn(start1, fixed = fixed, sumObs = TRUE,
      returnHessian = returnHessian, ...)
   if(slot(control, "printLevel") > 2) {
      cat("Initial function value:", f1, "\n")
   }
   if(any(is.na( f1))) {
      result <- list(code=100, message=maximMessage("100"),
                     iterations=0,
                     type=maximType)
      class(result) <- "maxim"
      return(result)
   }
   if(any(is.infinite( f1)) && sum(f1) > 0) {
                                        # we stop at +Inf but not at -Inf
      result <- list(code=5, message=maximMessage("5"),
                     iterations=0,
                     type=maximType)
      class(result) <- "maxim"
      return(result)
   }
   if( isTRUE( attr( f1, "gradBoth" ) ) ) {
      warning( "the gradient is provided both as attribute 'gradient' and",
         " as argument 'grad': ignoring argument 'grad'" )
   }
   if( isTRUE( attr( f1, "hessBoth" ) ) ) {
      warning( "the Hessian is provided both as attribute 'hessian' and",
         " as argument 'hess': ignoring argument 'hess'" )
   }
   G1 <- attr( f1, "gradient" )
   if(slot(control, "printLevel") > 2) {
      cat("Initial gradient value:\n")
      print(G1)
   }
   if(any(is.na(G1[!fixed]))) {
      stop("NA in the initial gradient")
   }
   if(any(is.infinite(G1[!fixed]))) {
      stop("Infinite initial gradient")
   }
   if(length(G1) != nParam) {
      stop( "length of gradient (", length(G1),
         ") not equal to the no. of parameters (", nParam, ")" )
   }
   H1 <- attr( f1, "hessian" )
   if(slot(control, "printLevel") > 3) {
      cat("Initial Hessian value:\n")
      print(H1)
   }
   if(length(H1) == 1) {
                           # Allow the user program to return a
                           # single NA in case of out of support or
                           # other problems
      if(is.na(H1))
          stop("NA in the initial Hessian")
   }
   if(any(is.na(H1[!fixed, !fixed]))) {
      stop("NA in the initial Hessian")
   }
   if(any(is.infinite(H1))) {
      stop("Infinite initial Hessian")
   }
   if( slot(control, "printLevel") > 1) {
      cat( "----- Initial parameters: -----\n")
      cat( "fcn value:",
      as.vector(f1), "\n")
      a <- cbind(start, G1, as.integer(!fixed))
      dimnames(a) <- list(nimed, c("parameter", "initial gradient",
                                          "free"))
      print(a)
      cat( "Condition number of the (active) hessian:",
          kappa( H1[!fixed, !fixed]), "\n")
      if( slot(control, "printLevel") > 3) {
         print( H1)
      }
   }
   lambda1 <- slot(control, "marquardt_lambda0")
   step <- 1
   ## ---------------- Main interation loop ------------------------
   repeat {
      if( iter >= slot(control, "iterlim")) {
         code <- 4; break
      }
      iter <- iter + 1
      if(!marquardt) {
         lambda1 <- 0
                           # assume the function is concave at start0
      }
      start0 <- start1
      f0 <- f1
      G0 <- G1
      if(any(is.na(G0[!fixed]))) {
         stop("NA in gradient (at the iteration start)")
      }
      H0 <- H1
      if(any(is.na(H0[!fixed, !fixed]))) {
         stop("NA in Hessian (at the iteration start)")
      }
      if(marquardt) {
         lambda1 <- lambda1/slot(control, "marquardt_lambdaStep")
                           # initially we try smaller lambda
                           # lambda1: current lambda for calculations
         H <- H0 - lambda1*I
      }
      else {
         step <- 1
         H <- H0
      }
      ## check whether hessian is positive definite
      aCount <- 0
                           # avoid inifinite number of attempts because of
                           # numerical problems
      while((me <-
         max.eigen( H[!fixed,!fixed,drop=FALSE])) >= -slot(control, "lambdatol") |
         (qRank <- qr(H[!fixed,!fixed], tol=slot(control, "qrtol"))$rank) < sum(!fixed)) {
                           # maximum eigenvalue -> negative definite
                           # qr()$rank -> singularity
         if(marquardt) {
            lambda1 <- lambda1*slot(control, "marquardt_lambdaStep")
         }
         else {
            lambda1 <- abs(me) + slot(control, "lambdatol") +
                min(abs(diag(H)[!fixed]))/1e7
                           # The third term corrects numeric singularity.  If diag(H) only contains large values,
                           # (H - (a small number)*I) == H because of finite precision
         }
         H <- (H - lambda1*I)
                           # could we multiply it with something like (for stephalving)
                           #     *abs(me)*lambdatol
                           # -lambda*I makes the Hessian (barely)
                           # negative definite.
                           # *me*lambdatol keeps the scale roughly
                           # the same as it was before -lambda*I
         aCount <- aCount + 1
         if(aCount > 100) {
                           # should be enough even in the worst case
            break
         }
      }
      amount <- vector("numeric", nParam)
      inv <- try(qr.solve(H[!fixed,!fixed,drop=FALSE],
                          G0[!fixed], tol=slot(control, "qrtol")))
      if(inherits(inv, "try-error")) {
                           # could not get the Hessian to negative definite
         samm <- list(theta0=start0, f0=f0, climb=amount)
         code <- 3
         break
      }
      amount[!fixed] <- inv
      start1 <- start0 - step*amount
                           # note: step is always 1 for Marquardt method
      f1 <- fn(start1, fixed = fixed, sumObs = TRUE,
         returnHessian = returnHessian, ...)
                           # The call calculates new function,
                           # gradient, and Hessian values
      ## Are we requested to fix some of the parameters?
      constPar <- attr(f1, "constPar")
      if(!is.null(constPar)) {
         if(any(is.na(constPar))) {
            stop("NA in the list of constants")
         }
         fixed <- rep(FALSE, nParam)
         fixed[constPar] <- TRUE
      }
      ## Are we asked to write in a new value for some of the parameters?
      if(is.null(newVal <- attr(f1, "newVal"))) {
         ## no ...
         if(marquardt) {
            stepOK <- lambda1 <= slot(control, "marquardt_maxLambda")
         }
         else {
            stepOK <- step >= slot(control, "steptol")
         }
         while( any(is.na(f1)) || ( ( sum(f1) < sum(f0) ) & stepOK)) {
                                        # We end up in a NA or a higher value.
                                        # try smaller step
            if(marquardt) {
               lambda1 <- lambda1*slot(control, "marquardt_lambdaStep")
               H <- (H0 - lambda1*I)
               amount[!fixed] <- qr.solve(H[!fixed,!fixed,drop=FALSE],
                                          G0[!fixed], tol=slot(control, "qrtol"))
            }
            else {
               step <- step/2
            }
            start1 <- start0 - step*amount
            if(slot(control, "printLevel") > 2) {
               if(slot(control, "printLevel") > 3) {
                  cat("Try new parameters:\n")
                  print(start1)
               }
               cat("function value difference", f1 - f0)
               if(marquardt) {
                  cat(" -> lambda", lambda1, "\n")
               }
               else {
                  cat(" -> step", step, "\n")
               }
            }
            f1 <- fn(start1, fixed = fixed, sumObs = TRUE,
               returnHessian = returnHessian, ...)
                           # WTF does the 'returnHessian' do here ?
            ## Find out the constant parameters -- these may be other than
            ## with full step
            constPar <- attr(f1, "constPar")
            if(!is.null(constPar)) {
               if(any(is.na(constPar))) {
                  stop("NA in the list of constants")
               }
               fixed[constPar] <- TRUE
               ## Any new values requested?
               if(!is.null(newVal <- attr(f1, "newVal"))) {
                  ## Yes.  Write them to parameters and go for
                  ## next iteration
                  start1[newVal$index] <- newVal$val
                  break;
               }
            }
         }
         if(marquardt) {
            stepOK <- lambda1 <= slot(control, "marquardt_maxLambda")
         }
         else {
            stepOK <- step >= slot(control, "steptol")
         }
         if(!stepOK) {
            # we did not find a better place to go...
            start1 <- start0
            f1 <- f0
            samm <- list(theta0=start0, f0=f0, climb=amount)
         }
      } else {
         ## Yes, indeed.  New values given to some of the params.
         ## Note, this may result in a lower function value,
         ## hence we do not check f1 > f0
         start1[newVal$index] <- newVal$val
         if( slot(control, "printLevel") > 0 ) {
            cat( "Keeping parameter(s) ",
               paste( newVal$index, collapse = ", " ),
               " at the fixed values ",
               paste( newVal$val, collapse = ", " ),
               ", as the log-likelihood function",
               " returned attributes 'constPar' and 'newVal'\n",
               sep = "" )
         }
      }
      G1 <- attr( f1, "gradient" )
      if(any(is.na(G1[!fixed]))) {
         cat("Iteration", iter, "\n")
         cat("Parameter:\n")
         print(start1)
         print(head(G1, n=30))
         stop("NA in gradient")
      }
      if(any(is.infinite(G1))) {
         code <- 6; break;
      }
      H1 <- attr( f1, "hessian" )
      if( slot(control, "printLevel") > 1) {
        cat( "-----Iteration", iter, "-----\n")
      }
      if(any(is.infinite(H1))) {
         code <- 7; break
      }
      if(slot(control, "printLevel") > 2) {
         cat( "lambda ", lambda1, " step", step, " fcn value:",
            formatC(as.vector(f1), digits=8, format="f"),  "\n")
         a <- cbind(amount, start1, G1, as.integer(!fixed))
         dimnames(a) <- list(names(start0), c("amount", "new param",
                                             "new gradient", "active"))
         print(a)
         if( slot(control, "printLevel") > 3) {
            cat("Hessian\n")
            print( H1)
         }
         if(!any(is.na(H1[!fixed, !fixed]))) {
            cat( "Condition number of the hessian:",
                kappa(H1[!fixed,!fixed,drop=FALSE]), "\n")
         }
      }
      if( step < slot(control, "steptol")) {
                           # wrong guess in step halving
         code <- 3; break
      }
      if(lambda1 > slot(control, "marquardt_maxLambda")) {
                           # wrong guess in Marquardt method
         code <- 3; break
      }
      if( sqrt( crossprod( G1[!fixed] ) ) < slot(control, "gradtol") ) {
         code <-1; break
      }
      if(is.null(newVal) && sum(f1) - sum(f0) < slot(control, "tol")) {
         code <- 2; break
      }
      if(is.null(newVal) && abs(sum(f1) - sum(f0)) <
         abs(slot(control, "reltol")*( sum(f1) + slot(control, "reltol")))) {
         code <- 2; break
      }
      if(any(is.infinite(f1)) && sum(f1) > 0) {
         code <- 5; break
      }
   }
   if( slot(control, "printLevel") > 0) {
      cat( "--------------\n")
      cat( maximMessage( code), "\n")
      cat( iter, " iterations\n")
      cat( "estimate:", start1, "\n")
      cat( "Function value:", f1, "\n")
   }
   names(start1) <- nimed
   F1 <- fn( start1, fixed = fixed, sumObs = FALSE,
      returnHessian = ( finalHessian == TRUE ), ... )
   G1 <- attr( F1, "gradient" )
   if(observationGradient(G1, length(start1))) {
      gradientObs <- G1
      colnames( gradientObs ) <- nimed 
      G1 <- colSums(as.matrix(G1 ))
   }
   else {
      gradientObs <- NULL
   }
   names( G1 ) <- nimed
   ## calculate (final) Hessian
   if(tolower(finalHessian) == "bhhh") {
      if(!is.null(gradientObs)) {
         hessian <- - crossprod( gradientObs )
         attr(hessian, "type") <- "BHHH"
      } else {
         hessian <- NULL
         warning("For computing the final Hessian by 'BHHH' method, the log-likelihood or gradient must be supplied by observations")
      }
   } else if( finalHessian != FALSE ) {
      hessian <- attr( F1, "hessian" )
   } else {
       hessian <- NULL
   }
   if( !is.null( hessian ) ) {
      rownames( hessian ) <- colnames( hessian ) <- nimed
   }

   ## remove attributes from final value of objective (likelihood) function
   attributes( f1 )$gradient <- NULL
   attributes( f1 )$hessian <- NULL
   attributes( f1 )$gradBoth <- NULL
   attributes( f1 )$hessBoth <- NULL
   ##
   result <-list(
       maximum = unname( drop( f1 ) ),
                  estimate=start1,
                  gradient=drop(G1),
                 hessian=hessian,
                  code=code,
       message=maximMessage( code),
                  last.step=samm,
                           # only when could not find a
                           # lower point
                  fixed=fixed,
       iterations=iter,
                  type=maximType)
   if( exists( "gradientObs" ) ) {
      result$gradientObs <- gradientObs
   }
   result <- c(result, control=control)
                           # attach the control parameters
   ##
   class(result) <- c("maxim", class(result))
   invisible(result)
}

returnCode.maxim <- function(x, ...)
    x$code
