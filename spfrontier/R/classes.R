# Copyright 2014 by Dmitry Pavlyuk <Dmitry.V.Pavlyuk@gmail.com>

#
# Definitions of S4 classes and methods
#


#' @title Model Estimation Results
#'
#' @description
#' \code{ModelEstimates} stores information about MLE estimates of a spatial stochastic frontier model
#' 
#' @details
#' \code{ModelEstimates} stores all parameter estimates and additional statistics, available after estimation of a spatial stochastic frontier model.
#' 
#' @slot coefficients estimated values of model parameters
#' @slot resultParams raw estimated values
#' @slot status model estimation status:\cr
#' 0 - Success\cr
#' 1 - Failed; convergence is not archieved\cr
#' 1000 - Failed; unexpected exception\cr
#' 1001 - Failed; Initial values for MLe can not be estimated\cr
#' 1002 - Failed; Maximum likelihood function is infinite\cr
#' 
#' @slot logL value of the log-likelihood function
#' @slot logLcalls information abour a number of log-likelihood function and its gradient function calls
#' @slot hessian Hessian matrix for estimated coefficients
#' @slot stdErrors standard errors of estimated coefficients
#' @slot residuals model residuals
#' @slot fitted model fitted values
#' @slot efficiencies estimates of efficiency values for sample observations
#' @name ModelEstimates-class
#' @rdname ModelEstimates-class
#' 

setClass("ModelEstimates", 
         representation(
             coefficients = "vector", 
             resultParams = "vector", 
             status = "numeric", 
             logL = "numeric", 
             logLcalls = "vector", 
             hessian = "matrix", 
             stdErrors = "vector", 
             residuals = "matrix", 
             fitted = "matrix", 
             efficiencies = "matrix")
)


#' Method \code{status} returns estimation status
#' @name ModelEstimates-class
#' @rdname ModelEstimates-class
#' @param object an object of ModelEstimates class
setGeneric("status",function(object){standardGeneric ("status")})

#' Method \code{resultParams} returns raw estimated coefficients
#' @name ModelEstimates-class
#' @rdname ModelEstimates-class
setGeneric("resultParams",function(object){standardGeneric ("resultParams")})

setGeneric("coefficients<-",function(object,value){standardGeneric("coefficients<-")})


setGeneric("residuals<-",function(object,value){standardGeneric("residuals<-")})


#' Method \code{hessian} returns Hessian matrix for estimated coefficients
#' @name ModelEstimates-class
#' @rdname ModelEstimates-class
setGeneric("hessian",function(object){standardGeneric("hessian")})
setGeneric("hessian<-",function(object,value){standardGeneric("hessian<-")})

#' Method \code{stdErrors} returns standard errors of estimated coefficients
#' @name ModelEstimates-class
#' @rdname ModelEstimates-class
setGeneric("stdErrors",function(object){standardGeneric("stdErrors")})
setGeneric("stdErrors<-",function(object,value){standardGeneric("stdErrors<-")})


setGeneric("fitted<-",function(object,value){standardGeneric("fitted<-")})

#' Method \code{efficiencies} returns efficiency estimates
#' @name ModelEstimates-class
#' @rdname ModelEstimates-class
setGeneric("efficiencies",function(object){standardGeneric ("efficiencies")})
setGeneric("efficiencies<-",function(object,value){standardGeneric("efficiencies<-")})



#' Method \code{show} prints estimated coefficients
#' @name ModelEstimates-class
#' @rdname ModelEstimates-class
#' @aliases show,ModelEstimates-method
setMethod("show", signature=signature(object='ModelEstimates'),
          function(object){
              cat("ModelEstimates\n", sep="")
              cat("Coefficients:\n")
              print(coefficients(object))
          }
)



#' Method \code{coefficients} returns estimated coefficients 
#' @rdname ModelEstimates-class
#' @aliases coefficients
setMethod("coefficients", signature=signature(object='ModelEstimates'),
          function(object){
              return(object@coefficients)
          }
)

#' @rdname ModelEstimates-class
#' @aliases resultParams
setMethod("resultParams", signature=signature(object='ModelEstimates'),
          function(object){
              return(unlist(object@resultParams))
          }
)

#' Method \code{fitted} returns model fitted values
#' @rdname ModelEstimates-class
#' @aliases fitted
setMethod("fitted", signature=signature(object='ModelEstimates'),
          function(object){
              return(object@fitted)
          }
)

#' @rdname ModelEstimates-class
#' @aliases efficiencies
setMethod("efficiencies", signature=signature(object='ModelEstimates'),
          function(object){
              return(object@efficiencies)
          }
)

#' Method \code{residuals} returns residuals
#' @rdname ModelEstimates-class
#' @aliases residuals
setMethod("residuals", signature=signature(object='ModelEstimates'),
          function(object){
              return(object@residuals)
          }
)

#' @rdname ModelEstimates-class
#' @aliases stdErrors
setMethod("stdErrors",signature=signature(object='ModelEstimates'),
          function(object) {
              return(object@stdErrors)
          }
)

#' @rdname ModelEstimates-class
#' @aliases hessian
setMethod("hessian",signature=signature(object='ModelEstimates'),
           function(object) {
               return(object@hessian)
           }
)


#' @rdname ModelEstimates-class
#' @aliases status
setMethod("status",signature=signature(object='ModelEstimates'), 
          function(object) {
              return(object@status)
        }
)

# Setters for the ModelEstimates class

setReplaceMethod("coefficients",signature=signature(object='ModelEstimates',value='vector'),
                 function(object, value){
                     object@coefficients <- value;
                     return(object)
                }
)


setReplaceMethod("fitted",signature=signature(object='ModelEstimates',value='matrix'), 
                 function(object, value){
                     object@fitted <- value;
                     return(object)
                }
)

setReplaceMethod("efficiencies",signature=signature(object='ModelEstimates',value='matrix'), 
                 function(object, value){
                     object@efficiencies <- value;
                     return(object)
                }
)

setReplaceMethod("residuals",signature=signature(object='ModelEstimates',value='matrix'), 
                 function(object, value){
                     object@residuals <- value;
                     return(object)
                }
)

setReplaceMethod("hessian",signature=signature(object='ModelEstimates',value='matrix'), 
                 function(object, value){
                     object@hessian <- value;
                     tryCatch({
                         OI <- solve(object@hessian)
                         v <- diag(OI)
                         se <- c()
                         for (i in v){
                             if (i>0){
                                 se <- c(se, sqrt(i))
                             } else{
                                 se <- c(se, NA)
                             }
                         }
                         object@stdErrors <- se
                     }, error = function(e){
                         logging(e$message, level="warn")
                     })
                     return(object)
                }
)

setReplaceMethod("stdErrors",signature=signature(object='ModelEstimates',value='vector'),
                 function(object, value){
                    object@stdErrors <- value;
                    return(object)
                }
)

# Summary of the estimated model

#' Method \code{summary} prints summary of the estimated model
#' @name ModelEstimates-class
#' @rdname ModelEstimates-class
#' @aliases summary,ModelEstimates-method
setMethod("summary", signature=signature(object='ModelEstimates'),
          function(object){
              cat("ModelEstimates\n", sep="")
              cat("Status: ", object@status,"\n")
              if (object@status == 0){
                  cat("Estimates:\n")
                  coef <- coefficients(object)
                  coef <- c(coef$beta,coef$rhoY, coef$sigmaV, coef$sigmaU, coef$rhoV, coef$rhoU, coef$mu)
                  zval <- coef/object@stdErrors
                  pval <- 2*(1-pnorm(abs(zval)))
                  digits <- max(5, .Options$digits - 3)
                  e <- cbind(format(signif(coef,digits=digits)), 
                             format(signif(object@stdErrors,digits=digits)), 
                             format(signif(zval,digits=digits)),
                             format.pval(pval,digits=digits),
                             sapply(pval,pvalMark))
                  rownames(e) <- names(coef)
                  colnames(e) <- c("Estimate","Std. Error", "z value", "Pr(>|z|)","")
                  print(e,quote = FALSE)
                  cat("---\n")
                  cat("Signif. codes:    0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1\n")
                  cat("Log-likelihood function value: ", object@logL,"\n")
                  cat("Log-likelihood function/gradient calls: ", unlist(object@logLcalls),"\n")
              }
          }
)
