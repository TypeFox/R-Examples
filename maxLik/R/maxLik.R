maxLik <- function(logLik, grad=NULL, hess=NULL, start,
                   method,
                   constraints=NULL,
                   ...) {
   ## Maximum Likelihood estimation.
   ##
   ## Newton-Raphson maximisation
   ## Parameters:
   ## logLik     log-likelihood function.  First argument must be the vector of parameters.
   ## grad       gradient of log-likelihood.  If NULL, numeric gradient is used.  Must return either
   ##               * vector, length=nParam
   ##               * matrix, dim=c(nObs, 1).  Treated as vector
   ##               * matrix, dim=c(nObs, nParam).  In this case the rows are simply
   ##                 summed (useful for maxBHHH).
   ## hess       Hessian function (numeric used if NULL)
   ## start      initial vector of parameters (eventually w/names)
   ## method     maximisation method (Newton-Raphson)
   ## constraints  constrained optimization: a list (see below)
   ## ...        additional arguments for the maximisation routine
   ##
   ## RESULTS:
   ## list of class c("maxLik", "maxim").  This is in fact equal to class "maxim", just the
   ## methods are different.
   ## maximum     function value at maximum
   ## estimate    the parameter value at maximum
   ## gradient        gradient
   ## hessian         Hessian
   ## code        integer code of success, depends on the optimization
   ##             method 
   ## message     character message describing the code
   ## type        character, type of optimization
   ##
   ##             there may be more components, depending on the choice of
   ##             the algorith.
   ##             
   argNames <-  c( "logLik", "grad", "hess", "start", "method",
                  "constraints" )
   checkFuncArgs( logLik, argNames, "logLik", "maxLik" )
   if( !is.null( grad ) ) {
      checkFuncArgs( grad, argNames, "grad", "maxLik" )
   }
   if( !is.null( hess ) ) {
      checkFuncArgs( hess, argNames, "hess", "maxLik" )
   }
   ## Constrained optimization.  We can two possibilities:
   ## * linear equality constraints
   ## * linear inequality constraints
   ##
   if(missing(method)) {
      if(is.null(constraints)) {
         method <- "nr"
      }
      else if(identical(names(constraints), c("ineqA", "ineqB"))) {
         if(is.null(grad))
             method <- "Nelder-Mead"
         else
             method <- "BFGS"
      }
      else
          method <- "nr"
   }
   maxRoutine <- switch(tolower(method),
                        "newton-raphson" =,
                        "nr" = maxNR,
                        "bfgs" = maxBFGS,
                        "bfgsr" =, "bfgs-r" = maxBFGSR,
                        "bhhh" = maxBHHH,
                        "conjugate-gradient" =,
                        "cg" = maxCG,
                        "nelder-mead" =,
                        "nm" = maxNM,
                        "sann" = maxSANN,
                        stop( "Maxlik: unknown maximisation method ", method )
                        )
   result <- maxRoutine(fn=logLik, grad=grad, hess=hess, start=start,
                        constraints=constraints,
                        ...)
   class(result) <- c("maxLik", class(result))
   result
}
