print.summary.maxLik <- function( x,
      digits = max( 3L, getOption("digits") - 3L ), ... ) {
   
   cat("--------------------------------------------\n")
   cat("Maximum Likelihood estimation\n")
   cat(maximType(x), ", ", nIter(x), " iterations\n", sep="")
   cat("Return code ", returnCode(x), ": ", returnMessage(x), "\n", sep="")
   if(!is.null(x$estimate)) {
      cat("Log-Likelihood:", x$loglik, "\n")
      cat(x$NActivePar, " free parameters\n")
      cat("Estimates:\n")
      printCoefmat( x$estimate, digits = digits )
   }
   if(!is.null(x$constraints)) {
      cat("\nWarning: constrained likelihood estimation.",
          "Inference is probably wrong\n")
      cat("Constrained optimization based on", x$constraints$type,
          "\n")
      if(!is.null(x$constraints$code))
         cat("Return code:", x$constraints$code, "\n")
                           # note: this is missing for 'constrOptim'
      if(!is.null(x$constraints$message))
         cat(x$constraints$message, "\n")
                           # note: this is missing for 'constrOptim'
      cat(x$constraints$outer.iterations,
          " outer iterations, barrier value",
          x$constraints$barrier.value, "\n")
   }
   cat("--------------------------------------------\n")
}

summary.maxLik <- function(object, eigentol=1e-12,... ) {
   ## object      object of class "maxLik"
   ## 
   ## RESULTS:
   ## list of class "summary.maxLik" with following components:
   ## maximum    : function value at optimum
   ## estimate   : estimated parameter values at optimum
   ## gradient   :           gradient at optimum
   ## code       : code of convergence
   ## message    : message, description of the code
   ## iterations : number of iterations
   ## type       : type of optimisation
   ##
   if(!inherits(object, "maxLik"))
       stop("'summary.maxLik' called on a non-'maxLik' object")
   ## Here we should actually coerce the object to a 'maxLik' object, dropping all the subclasses...
   ## Instead, we force the program to use maxLik-related methods
   result <- object$maxim
   nParam <- length(coef.maxLik(object))
   activePar <- activePar( object )
   if((object$code < 100) & !is.null(coef.maxLik(object))) {
                           # in case of infinity at initial values, the coefs are not provided
       t <- coef.maxLik(object)/stdEr.maxLik(object, eigentol=eigentol)
       p <- 2*pnorm( -abs( t))
       t[!activePar(object)] <- NA
       p[!activePar(object)] <- NA
       results <- cbind("Estimate"=coef.maxLik(object),
                        "Std. error"=stdEr.maxLik(object, eigentol=eigentol),
                        "t value"=t, "Pr(> t)"=p)
   } else {
     results <- NULL
   }
   summary <- list(maximType=object$type,
                   iterations=object$iterations,
                   returnCode=object$code,
                   returnMessage=object$message,
                   loglik=object$maximum,
                   estimate=results,
                   fixed=!activePar,
                   NActivePar=sum(activePar),
                   constraints=object$constraints)
   class(summary) <- "summary.maxLik"
   summary
}
