
### shoud move checkMaxControl to a separate file but how to do it?

setClassUnion("functionOrNULL", c("function", "NULL"))

checkMaxControl <- function(object) {
   ## check validity of MaxControl objects
   if(!inherits(object, "MaxControl")) {
      stop("'MaxControl' object required.  Currently '",
           class(object), "'")
   }
   ##
   errors <- character(0)
   ## Check length of componenents
   for(s in slotNames(object)) {
      if(s == "sann_cand") {
         if(length(slot(object, s)) > 1) {
            errors <- c(errors,
                        paste("'", s, "' must be either 'NULL' or ",
                              "a function of length 1, not of length ",
                              length(slot(object, s)), sep=""))
         }
      }
      else if(length(slot(object, s)) != 1) {
         errors <- c(errors,
                     paste("'", s, "' must be of length 1, not ",
                           length(slot(object, s)), sep=""))
      }
   }
   ##
   if(slot(object, "tol") < 0) {
      errors <- c(errors, paste("'tol' must be non-negative, not ",
                                slot(object, "tol"), sep=""))
   }
   if(slot(object, "reltol") < 0) {
      errors <- c(errors, paste("'reltol' must be non-negative, not ",
                                slot(object, "reltol"), sep=""))
   }
   if(slot(object, "gradtol") < 0) {
      errors <- c(errors, paste("'gradtol' must be non-negative, not",
                                slot(object, "gradtol")))
   }
   if(slot(object, "steptol") < 0) {
      errors <- c(errors, paste("'steptol' must be non-negative, not",
                                slot(object, "steptol")))
   }
   if(slot(object, "lambdatol") < 0) {
      errors <- c(errors, paste("'lambdatol' must be non-negative, not",
                                slot(object, "lambdatol")))
   }
   if(!pmatch(slot(object, "qac"), c("stephalving", "marquardt"))) {
      errors <- c(errors, paste("'qac' must be 'stephalving' or 'marquadt', not",
                                slot(object, "qac")))
   }
   if(slot(object, "qrtol") < 0) {
      errors <- c(errors, paste("'qrtol' must be non-negative, not",
                                slot(object, "qrtol")))
   }
   if(slot(object, "marquardt_lambda0") < 0) {
      errors <- c(errors, paste("'lambda0' must be non-negative, not",
                                slot(object, "lambda0")))
   }
   if(slot(object, "marquardt_lambdaStep") <= 1) {
      errors <- c(errors, paste("'lambdaStep' must be > 1, not",
                                slot(object, "lambdaStep")))
   }
   if(slot(object, "marquardt_maxLambda") < 0) {
      errors <- c(errors, paste("'maxLambda' must be non-negative, not",
                                slot(object, "maxLambda")))
   }
   ## NM
   if(slot(object, "nm_alpha") < 0) {
      errors <- c(errors, paste("Nelder-Mead reflection factor 'alpha' ",
                                "must be non-negative, not", slot(object, "nm_alpha")))
   }
   if(slot(object, "nm_beta") < 0) {
      errors <- c(errors, paste("Nelder-Mead contraction factor 'beta' ",
                                "must be non-negative, not", slot(object, "nm_beta")))
   }
   if(slot(object, "nm_gamma") < 0) {
      errors <- c(errors, paste("Nelder-Mead expansion factor 'gamma' ",
                                "must be non-negative, not", slot(object, "nm_gamma")))
   }
   ## SANN
   if(!inherits(slot(object, "sann_cand"), c("function", "NULL"))) { #
      errors <- c(errors, paste("'SANN_cand' must be either NULL or a function, not",
                                slot(object, "SANN_cand")))
   }
   if(slot(object, "sann_tmax") < 1) {
      errors <- c(errors, paste("SANN number of calculations at each temperature ",
                                "'tmax' ",
                                "must be positive, not", slot(object, "sann_tmax")))
   }
   ##
   if(slot(object, "iterlim") < 0) {
      errors <- c(errors, paste("'iterlim' must be non-negative, not",
                                slot(object, "iterlim")))
   }
   if(length(errors) > 0)
      return(errors)
   return(TRUE)
}

### MaxControls contains all control parameters for max* family
setClass("MaxControl",
         slots=representation(
             tol="numeric",
             reltol="numeric",
             gradtol="numeric",
             steptol="numeric",
                           #
             lambdatol="numeric",
             qrtol="numeric",
             ## Qadratic Approximation Control
             qac="character",
         marquardt_lambda0="numeric",
         marquardt_lambdaStep="numeric",
         marquardt_maxLambda="numeric",
         ## Optim Nelder-Mead:
         nm_alpha="numeric",
         nm_beta="numeric",
         nm_gamma="numeric",
         ## SANN
         sann_cand="functionOrNULL",
         sann_temp="numeric",
         sann_tmax="integer",
         sann_randomSeed="integer",
         ##
             iterlim="integer",
             ##
             printLevel="integer"),
         ##
         prototype=prototype(
             tol=1e-8,
             reltol=sqrt(.Machine$double.eps),
             gradtol=1e-6,
             steptol=1e-10,
                           #
             lambdatol=1e-6,
                           #
             qac="stephalving",
             qrtol=1e-10,
         marquardt_lambda0=1e-2,
         marquardt_lambdaStep=2,
         marquardt_maxLambda=1e12,
         ## Optim Nelder-Mead
         nm_alpha=1,
         nm_beta=0.5,
         nm_gamma=2,
         ## SANN
         sann_cand=NULL,
         sann_temp=10,
         sann_tmax=10L,
         sann_randomSeed=123L,
         ##
         iterlim=150L,
         printLevel=0L),
         ##
         validity=checkMaxControl)
