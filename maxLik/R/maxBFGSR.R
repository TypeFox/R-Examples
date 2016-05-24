
maxBFGSR <- function(fn, grad=NULL, hess=NULL, start, 
                     constraints=NULL,
                     finalHessian=TRUE,
                     fixed=NULL,
                     activePar=NULL,
                     control=NULL,
                     ...) {
   ## Newton-Raphson maximization
   ## Parameters:
   ## fn          - the function to be minimized.  Returns either scalar or
   ##               vector value with possible attributes 
   ##               constPar and newVal
   ## grad        - gradient function (numeric used if missing).  Must return either
   ##               * vector, length=nParam
   ##               * matrix, dim=c(nObs, 1).  Treated as vector
   ##               * matrix, dim=c(M, nParam), where M is arbitrary.  In this case the
   ##                 rows are simply summed (useful for maxBHHH).
   ## hess        - hessian function (numeric used if missing)
   ## start       - initial parameter vector (eventually w/names)
   ## ...         - extra arguments for fn()
   ##          The maxControl structure:
   ##       The stopping criteria
   ##     tol         - maximum allowed absolute difference between sequential values
   ##     reltol      - maximum allowed reltive difference (stops if < reltol*(abs(fn) + reltol)
   ##     gradtol     - maximum allowed norm of gradient vector
   ##     steptol     - minimum step size
   ##     iterlim     - maximum # of iterations
   ##     finalHessian  include final Hessian?  As computing final hessian does not carry any extra penalty for NR method, this option is
   ##                  mostly for compatibility reasons with other maxXXX functions.
   ##                  TRUE/something else  include
   ##                  FALSE                do not include
   ## activePar   - an index vector -- which parameters are taken as
   ##               variable (free).  Other paramters are treated as
   ##               fixed constants
   ## fixed         index vector, which parameters to keep fixed
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
   ##             activePar - logical vector, which parameters are active (not constant)
   ## activePar   logical vector, which parameters were treated as free (resp fixed)
   ## iterations  number of iterations
   ## type        "Newton-Raphson maximization"
   ##
   ## ------------------------------
   ## Add parameters from ... to control
   if(!inherits(control, "MaxControl")) {
      mControl <- addControlList(maxControl(), control)
   }
   else {
      mControl <- control
   }
   mControl <- addControlList(mControl, list(...), check=FALSE)
   ##
   argNames <- c(c( "fn", "grad", "hess", "start",
                 "activePar", "fixed", "control"),
                 openParam(mControl))
   checkFuncArgs( fn, argNames, "fn", "maxBFGSR" )
   if( !is.null( grad ) ) {
      checkFuncArgs( grad, argNames, "grad", "maxBFGSR" )
   }
   if( !is.null( hess ) ) {
      checkFuncArgs( hess, argNames, "hess", "maxBFGSR" )
   }
   ## establish the active parameters.  Internally, we just use 'activePar'
   fixed <- prepareFixed( start = start, activePar = activePar,
      fixed = fixed )
   ## chop off the control args from ... and forward the new ...
   dddot <- list(...)
   dddot <- dddot[!(names(dddot) %in% openParam(mControl))]
   cl <- list(start=start,
              finalHessian=finalHessian,
              fixed=fixed,
              control=mControl)
   if(length(dddot) > 0) {
      cl <- c(cl, dddot)
   }
   if(is.null(constraints)) {
      cl <- c(quote(maxBFGSRCompute),
              fn=logLikAttr,
              fnOrig = fn, gradOrig = grad, hessOrig = hess,
              cl)
      result <- eval(as.call(cl))
   } else {
      if(identical(names(constraints), c("ineqA", "ineqB"))) {
         stop("Inequality constraints not implemented for maxBFGSR")
      } else if(identical(names(constraints), c("eqA", "eqB"))) {
                           # equality constraints: A %*% beta + B = 0
         cl <- c(quote(sumt),
                 fn=fn, grad=grad, hess=hess,
                 maxRoutine=maxBFGSR,
                 constraints=list(constraints),
                 cl)
         result <- eval(as.call(cl))
      } else {
         stop("maxBFGSR only supports the following constraints:\n",
              "constraints=list(ineqA, ineqB)\n",
              "\tfor A %*% beta + B >= 0 linear inequality constraints\n",
              "current constraints:",
              paste(names(constraints), collapse=" "))
      }
   }
   return( result )
}
