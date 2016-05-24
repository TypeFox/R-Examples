maxOptim <- function(fn, grad, hess,
                    start, method, fixed,
                    constraints,
                     finalHessian=TRUE,
                    parscale,
                     control=maxControl(),
                    ...) {
   ## Wrapper of optim-based optimization methods
   ##
   ## finalHessian:   how (and if) to calculate the final Hessian:
   ##            FALSE   not calculate
   ##            TRUE    use analytic/numeric Hessian
   ##            bhhh/BHHH  use information equality approach
   ##
   if( method == "Nelder-Mead" ) {
      maxMethod <- "maxNM"
   } else {
      maxMethod <- paste( "max", method, sep = "" )
   }
   ##
   ## Add parameters from ... to control
   if(!inherits(control, "MaxControl")) {
      stop("'control' must be a 'MaxControl' object, created by 'maxControl()'")
   }
   control <- addControlList(control, list(...), check=FALSE)
   ## Any forbidden arguments in fn?
   argNames <- c( "fn", "grad", "hess", "start", "print.level", "iterlim",
      "constraints", "tol", "reltol", "parscale", "alpha", "beta", "gamma",
                 "cand", "temp", "tmax" )
   checkFuncArgs( fn, argNames, "fn", maxMethod )
   if( !is.null( grad ) ) {
      checkFuncArgs( grad, argNames, "grad", maxMethod )
   }
   if( !is.null( hess ) ) {
      checkFuncArgs( hess, argNames, "hess", maxMethod )
   }

   ## check argument 'fixed'
   fixed <- prepareFixed( start = start, activePar = NULL, fixed = fixed )

   message <- function(c) {
      switch(as.character(c),
             "0" = "successful convergence",
             "1" = "iteration limit exceeded",
             "10" = "degeneracy in Nelder-Mead simplex",
             "51" = "warning from the 'L-BFGS-B' method; see the corresponding component 'message' for details",
             "52" = "error from the 'L-BFGS-B' method; see the corresponding component 'message' for details"
             )
   }

   ## initialize variables for saving gradients provided as attributes
   ## and the corresponding parameter values
   lastFuncGrad <- NULL
   lastFuncParam <- NULL
   ## chop off the control args from '...' and forward the new '...'
   dddot <- list(...)
   dddot <- dddot[!(names(dddot) %in% openParam(control))]
                           # unfortunately now you have to do
                           # do.call(function, args, dddot) instead of just calling
                           # func(args, ...)
   ## strip possible SUMT parameters and call the function thereafter
   environment( callWithoutSumt ) <- environment()
   maximType <- paste( method, "maximization" )
   parscale <- rep(parscale, length.out=length(start))
   oControl <- list(trace=max(slot(control, "printLevel"), 0),
                    REPORT=1,
                    fnscale=-1,
                    reltol=slot(control, "tol"),
                    maxit=slot(control, "iterlim"),
                    parscale=parscale[ !fixed ],
                    alpha=slot(control, "nm_alpha"),
                    beta=slot(control, "nm_beta"),
                    gamma=slot(control, "nm_gamma"),
                    temp=slot(control, "sann_temp"),
                    tmax=slot(control, "sann_tmax")
                    )
   oControl$reltol <- slot(control, "reltol")
   argList <- list(theta=start,
                   fName="logLikFunc",
                   fnOrig = fn,
                   gradOrig = grad,
                   hessOrig = hess)
   if(length(dddot) > 0) {
      argList <- c(argList, dddot)
   }
   f1 <- do.call(callWithoutSumt, argList)
   if(is.na( f1)) {
      result <- list(code=100, message=maximMessage("100"),
                     iterations=0,
                     type=maximType)
      class(result) <- "maxim"
      return(result)
   }
   if(slot(control, "printLevel") > 2) {
      cat("Initial function value:", f1, "\n")
   }
   hasGradAttr <- !is.null( attr( f1, "gradient" ) )
   if( hasGradAttr && !is.null( grad ) ) {
      grad <- NULL
      warning( "the gradient is provided both as attribute 'gradient' and",
              " as argument 'grad': ignoring argument 'grad'" )
   }
   hasHessAttr <- !is.null( attr( f1, "hessian" ) )
   if( hasHessAttr && !is.null( hess ) ) {
      hess <- NULL
      warning( "the Hessian is provided both as attribute 'hessian' and",
              " as argument 'hess': ignoring argument 'hess'" )
   }
   if( method == "BFGS" ) {
      argList <- list(theta=start,
                      fName="logLikGrad",
                      fnOrig = fn,
                      gradOrig = grad,
                      hessOrig = hess)
      if(length(dddot) > 0) {
         argList <- c(argList, dddot)
      }
      G1 <- do.call(callWithoutSumt, argList)
      if(slot(control, "printLevel") > 2) {
         cat("Initial gradient value:\n")
         print(G1)
      }
      if(any(is.na(G1))) {
         stop("NA in the initial gradient")
      }
      if(any(is.infinite(G1))) {
         stop("Infinite initial gradient")
      }
      if(length(G1) != length(start)) {
         stop( "length of gradient (", length(G1),
              ") not equal to the no. of parameters (", length(start), ")" )
      }
   }

   ## function to return the gradients (BFGS, CG) or the new candidate point (SANN)
   if( method == "BFGS" ) {
      gradOptim <- logLikGrad
   } else if( method == "SANN" ) {
      if( is.null(slot(control, "sann_cand") ) ) {
         gradOptim <- NULL
      } else {
         gradOptim <- function( theta, fnOrig, gradOrig, hessOrig,
               start, fixed, ... ) {
            return(control@sann_cand( theta, ... ) )
         }
      }
   } else if( method == "CG" ) {
      gradOptim <- logLikGrad
   } else if( method == "Nelder-Mead" ) {
      gradOptim <- NULL
   } else {
      stop( "internal error: unknown method '", method, "'" )
   }
   ## A note about return value:
   ## We can the return from 'optim' in a object of class 'maxim'.
   ## However, as 'sumt' already returns such an object, we return the
   ## result of 'sumt' directly, without the canning
   if(is.null(constraints)) {
      cl <- list(quote(optim),
                 par = start[ !fixed ],
                 fn = logLikFunc,
                 control = oControl,
                 method = method,
                 gr = gradOptim,
                 fnOrig = fn,
                 gradOrig = grad,
                 hessOrig = hess,
                 start = start,
                 fixed = fixed)
      if(length(dddot) > 0) {
         cl <- c(cl, dddot)
      }
      result <- eval(as.call(cl))
      resultConstraints <- NULL
   }
   else {
      ## linear equality and inequality constraints
                           # inequality constraints: A %*% beta + B >= 0
      if(identical(names(constraints), c("ineqA", "ineqB"))) {
         nra <- nrow(constraints$ineqA)
         nrb <- nrow(as.matrix(constraints$ineqB))
         ncb <- ncol(as.matrix(constraints$ineqB))
         if(ncb != 1) {
            stop("Inequality constraint B must be a vector ",
                 "(or Nx1 matrix).  Currently ", ncb, " columns")
         }
         if(length(dim(constraints$ineqA)) != 2) {
            stop("Inequality constraint A must be a matrix\n",
                 "Current dimension", dim(constraints$ineqA))
         }
         if(ncol(constraints$ineqA) != length(start)) {
            stop("Inequality constraint A must have the same ",
                 "number of columns as length of the parameter.\n",
                 "Currently ", ncol(constraints$ineqA),
                 " and ", length(start), ".")
         }
         if(ncol(constraints$ineqA) != length(start)) {
            stop("Inequality constraint A cannot be matrix multiplied",
                 " with the start value.\n",
                 "A is a ", nrow(constraints$ineqA), "x",
                 ncol(constraints$ineqA), " matrix,",
                 " start value has lenght ", length(start))
         }
         if(nra != nrb) {
            stop("Inequality constraints A and B suggest different number ",
                 "of constraints: ", nra, " and ", nrb)
         }
         cl <- list(quote(constrOptim2),
                    theta = start,
                    f = logLikFunc, grad = gradOptim,
                    ineqA=constraints$ineqA,
                    ineqB=constraints$ineqB,
                    control=oControl,
                    method = method, fnOrig = fn, gradOrig = grad,
                    hessOrig = hess, fixed = fixed, start=start)
                           # 'start' argument is needed for adding fixed parameters later in the call chain
         if(length(dddot) > 0) {
            cl <- c(cl, dddot)
         }
         result <- eval(as.call(cl))
         resultConstraints <- list(type="constrOptim",
                                   barrier.value=result$barrier.value,
                                   outer.iterations=result$outer.iterations
                                   )
      }
      else if(identical(names(constraints), c("eqA", "eqB"))) {
                           # equality constraints: A %*% beta + B = 0
         argList <- list(fn=fn, grad=grad, hess=hess,
                        start=start, fixed = fixed,
                        maxRoutine = get( maxMethod ),
                        constraints=constraints,
                         parscale = parscale,
                         control=control)
                           # recursive evaluation-> pass original (possibly
                           # supplemented) control
         if(length(dddot) > 0) {
            argList <- c(argList, dddot)
         }
         result <- do.call( sumt, argList[ !sapply( argList, is.null ) ] )
         return(result)
                           # this is already maxim object
      }
      else {
         stop( maxMethod, " only supports the following constraints:\n",
              "constraints=list(ineqA, ineqB)\n",
              "\tfor A %*% beta + B >= 0 linear inequality constraints\n",
              "current constraints:",
              paste(names(constraints), collapse=" "))
      }
   }

   # estimates (including fixed parameters)
   estimate <- start
   estimate[ !fixed ] <- result$par
   ## Calculate the final gradient
   argList <- list(estimate, "logLikGrad",
              fnOrig = fn, gradOrig = grad, hessOrig = hess, sumObs = FALSE)
   if(length(dddot) > 0) {
      argList <- c(argList, dddot)
   }
   gradient <- do.call(callWithoutSumt, argList)
   if(observationGradient(gradient, length(start))) {
      gradientObs <- gradient
      gradient <- colSums(as.matrix(gradient ))
   }
   else {
      gradientObs <- NULL
   }
   ## calculate (final) Hessian
   if(tolower(finalHessian) == "bhhh") {
      if(!is.null(gradientObs)) {
         hessian <- - crossprod( gradientObs )
         attr(hessian, "type") <- "BHHH"
      } else {
         hessian <- NULL
         warning("For computing the final Hessian by 'BHHH' method, the log-likelihood or gradient must be supplied by observations")
      }
   } else if(finalHessian != FALSE) {
      argList <- list( estimate, fnOrig = fn,  gradOrig = grad,
                            hessOrig = hess)
      if(length(dddot) > 0) {
         argList <- c(argList, dddot)
      }
      hessian <- as.matrix( do.call(logLikHess, argList) )
   } else {
       hessian <- NULL
   }
   if( !is.null( hessian ) ) {
      rownames( hessian ) <- colnames( hessian ) <- names( estimate )
   }

   result <- list(
                   maximum=result$value,
                   estimate=estimate,
                   gradient=drop(gradient),
                           # ensure the final (non-observation) gradient is just a vector
                   hessian=hessian,
                   code=result$convergence,
                   message=paste(message(result$convergence), result$message),
                   last.step=NULL,
                   fixed = fixed,
                   iterations=result$counts[1],
                   type=maximType,
                  constraints=resultConstraints
                  )
   if( exists( "gradientObs" ) ) {
      result$gradientObs <- gradientObs
   }
   result <- c(result, control=control)
                           # attach the control parameters
   class(result) <- "maxim"
   return(result)
}
