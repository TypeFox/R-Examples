print.summary.maxim <- function( x, ... ) {
   summary <- x
   cat("--------------------------------------------\n")
   cat(summary$type, "\n")
   cat("Number of iterations:", summary$iterations, "\n")
   cat("Return code:", summary$code, "\n")
   cat(summary$message, "\n")
   if(!is.null(summary$unsucc.step)) {
      cat("Last (unsuccessful) step: function value", summary$unsucc.step$value,
         "\n")
      print(summary$unsucc.step$parameters)
   }
   if(!is.null(summary$estimate)) {
      cat("Function value:", summary$maximum, "\n")
      cat("Estimates:\n")
      print(summary$estimate)
      if(!is.null(summary$hessian)) {
         cat("Hessian:\n")
         print(summary$hessian)
      }
   }
   if(!is.null(summary$constraints)) {
      cat("\nConstrained optimization based on", summary$constraints$type,
          "\n")
      if(!is.null(summary$constraints$code))
         cat("Return code:", summary$constraints$code, "\n")
                           # note: this is missing for 'constrOptim'
      if(!is.null(summary$constraints$message))
         cat(summary$constraints$message, "\n")
                           # note: this is missing for 'constrOptim'
      cat(summary$constraints$outer.iterations,
          " outer iterations, barrier value",
          summary$constraints$barrier.value, "\n")
   }
   cat("--------------------------------------------\n")
}

summary.maxim <- function(object, hessian=FALSE, unsucc.step=FALSE,
   ... ) {
   ## The object of class "maxim" should include following components:
   ## maximum    : function value at optimum
   ## estimate   : matrix, estimated parameter values and gradient at optimum
   ## hessian    :           hessian
   ## code       : code of convergence
   ## message    : message, description of the code
   ## last.step  : information about last step, if unsuccessful
   ## iterations : number of iterations
   ## type       : type of optimisation
   ##
   nParam <- length(object$estimate)
   if(object$code == 3 & unsucc.step) {
      a <- cbind(object$last.step$theta0, object$last.step$theta1)
      dimnames(a) <- list(parameter=object$names,
                        c("current par", "new par"))
      unsucc.step <- list(value=object$last.step$f0,
                        parameters=a)
   } else {
      unsucc.step <- NULL
   }
   estimate <- cbind("estimate"=object$estimate, "gradient"=object$gradient)
   if(hessian) {
      H <- object$hessian
   }
   else {
      H <- NULL
   }
   summary <- list(maximum=object$maximum,
                   type=object$type,
                  iterations=object$iterations,
                  code=object$code,
                  message=object$message,
                   unsucc.step=unsucc.step,
                  estimate=estimate,
                   hessian=H,
                   constraints=object$constraints)
   class(summary) <- c("summary.maxim", class(summary))
   summary
}
