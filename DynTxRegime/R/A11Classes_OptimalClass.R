#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                              Class OptimalClass                              #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# contains all relevant results of the doubly robust dynamic tx regime         #
# methodology.                                                                 #
#                                                                              #
#  value   : E(g^{opt} )                                                        #
#                                                                              #
#  classif : The classification scheme object.                                 #
#                                                                              #
#  propen  : The fit object for each propensity of tx model.                   #
#                                                                              #
#  outcome : The fit object for the conditional expectation model.             #
#                                                                              #
#  optTx   : vector of optimal treatments for training set                     #
#                                                                              #
#  call    : matching call                                                     #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
if(!isClass("QLearnEst or NULL")){
  setClassUnion("QLearnEst or NULL", 
                members = c("QLearnEst","NULL"))
}

setClass(Class = "OptimalClass",
         slots = c( estVal = "numeric",
                   classif = "modelObjFit",
                    propen = "PropenFit",
                   outcome = "QLearnEst or NULL",
                     optTx = "factor",
                      call = "call"),
         contains = "DynTxRegime")

setMethod(f = "Call", 
          signature = c(object = "OptimalClass"), 
          definition = function(object, ...){ 
                         return(object@call) 
                       } )

setMethod(f = "Classif",
          signature = c(object="OptimalClass"),
          definition = function(object, ...){ 
                         return( object@classif )  
                       } )

setMethod(f = "Coef",    
          signature = c(object = "OptimalClass"), 
          definition = function(object, ...){
                         result <- list()
                         result[[ "propen" ]] <- 
                           Coef(object = object@propen, ...)
                         if( is(object@outcome, "NULL") ) return( result )
                         result[[ "outcome" ]] <- 
                           Coef(object = object@outcome, ...)
                         return(result)
                       } )

setMethod(f = "Est",    
          signature = c(x = "OptimalClass"), 
          definition = function(x, ...){
                         return( x@estVal )
                       } )

setMethod(f = "FitObject",    
          signature = c(object = "OptimalClass"), 
          definition = function(object, ...){
                         result <- list()
                         result[[ "classif" ]] <- 
                           FitObject(object = object@classif)
                         result[[ "propen" ]] <- 
                           FitObject(object = object@propen, ...)
                         if( is(object@outcome, "NULL") ) return( result )
                         result[[ "outcome" ]] <- 
                           FitObject(object = object@outcome, ...)
                         return(result)
                       } )

setMethod(f = "FittedCont",    
          signature = c(object = "OptimalClass"), 
          definition = function(object, ...){
                         return( FittedCont(object = object@outcome, ...) )
                       } )

setMethod(f = "FittedMain",    
          signature = c(object = "OptimalClass"), 
          definition = function(object, ...){
                         return( FittedMain(object = object@outcome, ...) )
                       } )

setMethod(f = "Genetic",    
          signature = c(object = "OptimalClass"), 
          definition = function(object, ...){
                         return( NULL )
                       } )

setMethod(f = "ModelObjectFit",    
          signature = c(object = "OptimalClass"), 
          definition = function(object, ...){
                         result <- list()
                         result[[ "classif" ]] <- 
                           ModelObjectFit(object = object@classif, ...)
                         result[[ "propen" ]] <- 
                           ModelObjectFit(object = object@propen, ...)
                         if( is(object@outcome) ) return( result )
                         result[[ "outcome" ]] <- 
                           ModelObjectFit(object = object@outcome, ...)
                         return(result)
                       } )

setMethod(f = "OptTx",  
          signature = c(x = "OptimalClass",
                        newdata = "missing"), 
          definition = function(x, newdata, ...){ 
                         return( x@optTx ) 
                       } )

setMethod(f = "OptTx",
          signature = c(x = "OptimalClass", 
                        newdata = "data.frame"),
          definition = function (x, newdata,...){
                         pred1 <- predict(object = x@classif, newdata = newdata)
                         opt <- as.numeric(as.character(pred1))
                         return( opt )
                       } )

setMethod(f = "Outcome", 
          signature = c(object = "OptimalClass"), 
          definition = function(object, ...){
                         if( is(object@outcome, "NULL") ) return(NULL)
                         return( FitObject(object = object@outcome, ...) )
                       } )

setMethod(f = "Plot",
          signature = c(x = "OptimalClass"),
          definition = function(x, suppress=FALSE, ...){

                         argList <- list(...)
                         if( !suppress ) {
                           nms <- "Classification"
                           if( is(argList[[ "main" ]], "NULL") ) {
                             argList[[ "main" ]] <- nms
                           } else if( is(argList[[ "sub" ]], "NULL") ) {
                             argList[[ "sub" ]] <- nms
                           } else {
                             argList[[ "sub" ]] <- paste(argList[[ "sub" ]], 
                                                         " (",nms,")", sep="")
                           }
                         }
                         argList[[ "x" ]] <- x@classif
                         argList[[ "suppress" ]] <- suppress
                         do.call(what = Plot, args = argList)

                         argList <- list(...)
                         if( !suppress ) {
                           nms <- "Propensity for Treatment"
                           if( is(argList[[ "main" ]], "NULL") ) {
                             argList[[ "main" ]] <- nms
                           } else if( is(argList[[ "sub" ]], "NULL") ) {
                             argList[[ "sub" ]] <- nms
                           } else {
                             argList[[ "sub" ]] <- paste(argList[[ "sub" ]], 
                                                         " (",nms,")", sep="")
                           }
                         }
                         argList[[ "x" ]] <- x@propen
                         argList[[ "suppress" ]] <- suppress
                         do.call(what = Plot, args = argList)

                         if( is(x@outcome, "NULL") ) return()

                         argList <- list(...)
                         if( !suppress ) {
                           nms <- "Outcome Regression"
                           if( is(argList[[ "main" ]], "NULL") ) {
                             argList[[ "main" ]] <- nms
                           } else if( is(argList[[ "sub" ]], "NULL") ) {
                             argList[[ "sub" ]] <- nms
                           } else {
                             argList[[ "sub" ]] <- paste(argList[[ "sub" ]], 
                                                       " (",nms,")", sep="")
                           }
                         }
                         argList[[ "x" ]] <- x@outcome
                         argList[[ "suppress" ]] <- suppress
                         do.call(what = Plot, args = argList)
                       } )

setMethod(f = "Print",    
          signature = c(x = "OptimalClass"), 
          definition = function(x, ...){
                         cat("\nClassification\n")
                         Print(x@classif)
                         cat("\nPropensity for Treatment\n")
                         Print(x@propen, ...)
                         if( !is(x@outcome, "NULL") ) {
                           cat("\nOutcome Regression\n")
                           Print(x@outcome, ...)
                         }
                         cat("\nEstimated Value:", x@estVal,"\n\n")
                         return()
                       } )

setMethod(f = "Propen", 
          signature = c(object = "OptimalClass"), 
          definition = function(object, ...){
                         return( FitObject(object = object@propen, ...) )
                       } )

setMethod(f = "RegimeCoef",    
          signature = c(object = "OptimalClass"), 
          definition = function(object, ...){
                         return( NULL )
                       } )

setMethod(f = "Residuals",    
          signature = c(object = "OptimalClass"), 
          definition = function(object, ...){
                         result <- list()
                         result[[ "propen" ]] <- 
                           Residuals(object = object@propen, ...)
                         if( is(object@outcome, "NULL") ) return( result )
                         result[[ "outcome" ]] <- 
                           Residuals(object = object@outcome, ...)
                         return(result)
                       } )

setMethod(f = "Show",    
          signature = c(object = "OptimalClass"), 
          definition = function(object, ...){
                         cat("Classification\n")
                         show(object@classif)
                         cat("Propensity for Treatment\n")
                         Show(object@propen, ...)
                         if( !is(object@outcome, "NULL") ) {
                           cat("Outcome Regression\n")
                           Show(object@outcome, ...)
                         }
                         cat("Estimated Value:", object@estVal,"\n")
                         return()
                       } )

setMethod(f = "StdDev",    
          signature = c(object = "OptimalClass"), 
          definition = function(object, ...){
                         return( NULL )
                       } )

setMethod(f = "Step", 
          signature = c(object = "OptimalClass"), 
          definition = function(object, ...){
                         return("OptimalClass")
                       } )

setMethod(f = "MySummary",    
          signature = c(object = "OptimalClass"), 
          definition = function(object, ...){
                         result <- list()
                         result[[ "classif" ]] <- 
                           MySummary(object = object@classif, ...)
                         result[[ "propen" ]] <- 
                           MySummary(object = object@propen, ...)
                         if( is(object@outcome, "NULL") ) return( result )
                         result[[ "outcome" ]] <- 
                           MySummary(object = object@outcome, ...)
                         return(result)
                       } )


