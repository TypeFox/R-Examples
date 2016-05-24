openParam <- function(object) {
   ## Return character list of 'open parameters', parameters that can
   ## be supplied to max* outside of 'control' list
   ## 
   if(!inherits(object, "MaxControl")) {
      stop("'MaxControl' object required.  Currently ",
           class(object))
   }
   c("tol",
     "reltol",
     "gradtol",
     "steptol",
     #
     "lambdatol",
     ## Qadratic Approximation Control
     "qac",
     "qrtol",
     "lambda0",
     "lambdaStep",
     "maxLambda",
     ## optim Nelder-Mead
     "alpha", "beta", "gamma",
     ## SANN (open versions)
     "cand", "temp", "tmax", "random.seed",
     ##
     "iterlim",
     "printLevel", "print.level")
}
