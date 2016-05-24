#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                                 CLASS Regime                                 #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# contains all relevant information for a regime function                      #
#                                                                              #
#   nvars  : the number of parameters to be estimated in the regime.           #
#                                                                              #
#   vnames : the names of the parameters to be estimated in the regime         #
#                                                                              #
#   func   : the user specified function that defines the regime               #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

setClass("Regime",
         slots = c( nVars = "integer",
                   vNames = "character",
                     func = "function") )   

setClass(Class = "RegimeList",
         contains = "List" )

setMethod(f = "NVars", 
          signature = c(object="Regime"), 
          definition = function(object, ...){
                         return( object@nVars ) 
                       } )
setMethod(f = "VNames",  
          signature = c(object="Regime"), 
          definition = function(object, ...){ 
                         return( object@vNames ) 
                       } )
setMethod(f = "RegFunc",   
          signature = c(object="Regime"), 
          definition = function(object, ...){ 
                         return( object@func ) 
                       } )

if(!isClass("Regime or RegimeList")){
  setClassUnion("Regime or RegimeList", 
                members = c("Regime", "RegimeList") )
}

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                               Class OptimalSeq                               #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# contains all relevant results of the doubly robust dynamic tx regime         #
# methodology.                                                                 #
#                                                                              #
#  est       : Parameter estimates for the optimal tx regime.                  #
#                                                                              #
#  value     : Value of the AIPWE estimator at the estimated optimal tx regime #
#                                                                              #
#  gaFit     : Object returned by genoud. See ?rgenoud # for details.          #
#                                                                              #
#  propen : Fit object for each propensity of tx model.                        #
#                                                                              #
#  outcome  : Fit object for each conditional expectation model.               #
#              If refit=FALSE, these are the fitted Q-functions.               #
#              If refit=TRUE, these are the parameter estimates based on the   #
#              estimated optimal tx regime.                                    #
#                                                                              #
#  regimes : regime information provided by user                               #
#                                                                              #
#  optTx : matrix of optimal treatments at each dp for the training set        #
#                                                                              #
#  call      : matching call                                                   #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
if(!isClass("QLearnEst or QLearnEstList or NULL")){
  setClassUnion("QLearnEst or QLearnEstList or NULL", 
                members = c("QLearnEst","QLearnEstList","NULL"))
}

setClass("OptimalSeq",
         slots = c( varEst = "list",
                    estVal = "numeric",
                   genetic = "list",
                    propen = "PropenFit or PropenFitList",
                   outcome = "QLearnEst or QLearnEstList or NULL",
                   regimes = "Regime or RegimeList",
                     optTx = "matrix",
                    txInfo = "TxInfo or TxInfoList",
                      call = "call"),
         contains = "DynTxRegime")

setMethod(f = "Call",    
          signature = c(object = "OptimalSeq"), 
          definition = function(object, ...){  
                         return( object@call )  
                       } )

setMethod(f = "Coef",    
          signature = c(object = "OptimalSeq"), 
          definition = function(object, ...){
                         result <- list()
                         result[[ "regimes" ]] <- object@varEst
                         result[[ "propen" ]] <- 
                           Coef(object = object@propen, ...)
                         if( is(object@outcome, "NULL") ) return(result)
                         result[[ "outcome" ]] <- 
                           Coef(object = object@outcome, ...)
                         return(result)
                       } )

setMethod(f = "Est",    
          signature = c(x = "OptimalSeq"), 
          definition = function(x, ...){
                         return( x@estVal )
                       } )

setMethod(f = "FitObject",    
          signature = c(object = "OptimalSeq"), 
          definition = function(object, ...){
                         result <- list()
                         result[[ "genetic" ]] <- object@genetic
                         result[[ "propen" ]] <- 
                           FitObject(object = object@propen, ...)
                         if( is(object@outcome, "NULL") ) return(result)
                         result[[ "outcome" ]] <- 
                           FitObject(object = object@outcome, ...)
                         return(result)
                       } )

setMethod(f = "FittedCont",    
          signature = c(object = "OptimalSeq"), 
          definition = function(object, ...){
                         return( FittedCont(object = object@outcome) )
                       } )

setMethod(f = "FittedMain",    
          signature = c(object = "OptimalSeq"), 
          definition = function(object, ...){
                         return( FittedMain(object = object@outcome) )
                       } )

setMethod(f = "Genetic",    
          signature = c(object = "OptimalSeq"), 
          definition = function(object, ...){  
                         return( object@genetic ) 
                       } )

setMethod(f = "ModelObjectFit",    
          signature = c(object = "OptimalSeq"), 
          definition = function(object, ...){
                         result <- list()
                         result[[ "propen" ]] <- 
                           ModelObjectFit(object = object@propen, ...)
                         if( is(object@outcome, "NULL") ) return( result )
                         result[[ "outcome" ]] <- 
                           ModelObjectFit(object = object@outcome, ...)
                         return(result)
                       } )

setMethod(f = "OptTx",    
          signature = c(x = "OptimalSeq",
                        newdata = "missing"), 
          definition = function(x, newdata, dp=1L, ...){  
                         return( x@optTx[,dp] ) 
                       } )

setMethod(f = "OptTx",
          signature = c(x = "OptimalSeq", 
                        newdata = "data.frame"),
          definition = function (x, newdata, dp=1L, ...){

                         regs <- x@regimes
                         pars <- x@varEst

                         optTx <- matrix(NA,
                                         nrow=nrow(newdata),
                                         ncol=1L,
                                         dimnames=list(NULL, 
                                                  paste("dp=", dp)))

                         if( is(regs, 'RegimeList') ) {
                           rr <- RegFunc(regs[[dp]])
                           pp <- pars[[dp]]
                           tx <- TxName(x@txInfo[[dp]])
                         } else if( is(regs, 'Regime') ) {
                           rr <- RegFunc(regs)
                           pp <- pars[[1L]]
                           tx <- TxName(x@txInfo)
                         } else {
                           DeveloperError("bad regime class", "OptTx")
                         }

                         argList <- list()
                         nms <- names(formals(rr))
                         np <- length(nms)
                         parNames <- nms[1L:{np-1L}]
                         dataName <- nms[np]
                         for(j in 1L:length(pp)) {
                           argList[[ parNames[j] ]] <- pp[j]
                         }

                         argList[[dataName]] <- newdata
                         reg.g <- do.call(what = rr, args = argList)

                         optTx[,1L] <- reg.g

                         return( optTx )
                       } )

setMethod(f = "Outcome",    
          signature = c(object = "OptimalSeq"), 
          definition = function(object, ...){
                         if( is(object@outcome, "NULL") ) return(NULL)
                         return( FitObject(object = object@outcome) )
                       } )

setMethod(f = "Plot",
          signature = c(x = "OptimalSeq"),
          definition = function(x, suppress=FALSE, ...){
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
          signature = c(x = "OptimalSeq"), 
          definition = function(x, ...){
                         cat("\nGenetic Algorithm\n")
                         show(x@genetic)
                         cat("\nPropensity for Treatment\n")
                         Show(x@propen, ...)
                         if( !is(x@outcome, "NULL") ) {
                           cat("\nOutcome Regression\n")
                           Show(x@outcome, ...)
                         }
                         cat("\nRegime Parameters:\n")
                         print(RegimeCoef(x))
                         cat("\nEstimated Value:", Est(x),"\n")
                         return()
                       } )

setMethod(f = "Propen",    
          signature = c(object = "OptimalSeq"), 
          definition = function(object, ...){
                         return( FitObject(object = object@propen) )
                       } )

setMethod(f = "RegimeCoef",    
          signature = c(object = "OptimalSeq"), 
          definition = function(object, ...){
                         me <- object@varEst
                         if( is(me, "list") ) {
                           if( length(me) == 1L ) me <- me[[1]]
                         }
                         return( me )
                       } )

setMethod(f = "Regimes",    
          signature = c(object = "OptimalSeq"), 
          definition = function(object, ...){ 
                         return( object@regimes ) 
                       } )

setMethod(f = "Residuals",    
          signature = c(object = "OptimalSeq"), 
          definition = function(object, ...){
                         result <- list()
                         result[[ "propen" ]] <- 
                           Residuals(object = object@propen, ...)
                         if( is(object@outcome, "NULL") ) return( result )
                         result[[ "outcome" ]] <- 
                           Residuals(object = object@outcome, ...)
                         return( result )
                       } )

setMethod(f = "Show",    
          signature = c(object = "OptimalSeq"), 
          definition = function(object, ...){
                         cat("\nGenetic Algorithm\n")
                         show(object@genetic)
                         cat("\nPropensity for Treatment\n")
                         Show(object@propen, ...)
                         if( !is(object@outcome, "NULL") ) {
                           cat("\nOutcome Regression\n")
                           Show(object@outcome, ...)
                         }
                         cat("\nRegime Parameters:\n")
                         print(RegimeCoef(object))
                         cat("\n\nEstimated Value:", Est(object),"\n")
                         return()
                       } )

setMethod(f = "Step", 
          signature = c(object = "OptimalSeq"), 
          definition = function(object, ...){ 
                         return("OptimalSeq") 
                       } )

setMethod(f = "MySummary",    
          signature = c(object = "OptimalSeq"), 
          definition = function(object, ...){
                         result <- list()
                         result[[ "propen" ]] <- 
                           MySummary(object = object@propen, ...)
                         if( is(object@outcome, "NULL") ) return( result )
                         result[[ "outcome" ]] <- 
                           MySummary(object = object@outcome, ...)
                         return(result)
                       } )

