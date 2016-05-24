#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                               CLASS QLearnEst                                #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# contains results obtained by fitting models for a stage of Q-learning        #
#                                                                              #
#   fitObj     : an object of class simpleFit, iterateFit, a list of objects of#
#                class simpleFit, or a list of objects of class iterateFit     #
#                                                                              #
#   qFunctions : matrix of Q-function values                                   #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setClass(Class = "QLearnEst",
         slots = c(   fitObj = "SimpleFit, IterateFit, SubsetFitList",
                   valueFunc = "numeric") )

setClass(Class = "QLearnEstList",
         contains = "dpList" )

setMethod(f = "[[<-", 
          signature = c(x = "QLearnEstList"), 
          definition = function(x,i,value){
                         x@loo[[i]] <- value
                         return(x)
                       } )

setMethod(f = "Base",    
          signature = c(object = "QLearnEst"), 
          definition = function(object, ...){ 
                         return( Base(object@fitObj) ) 
                       } )

setMethod(f = "Coef", 
          signature = c(object = "QLearnEst"), 
          definition = function(object, ...){
                         return( Coef(object = object@fitObj, ...) )
                       } )

setMethod(f = "FitObject", 
          signature = c(object = "QLearnEst"), 
          definition = function(object, ...){
                         return( FitObject(object = object@fitObj, ...) )
                       } )

setMethod(f = "FittedCont", 
          signature = c(object = "QLearnEst"), 
          definition = function(object, ...){
                         return( FittedCont(object = object@fitObj, ...) )
                       } )

setMethod(f = "FittedMain", 
          signature = c(object = "QLearnEst"), 
          definition = function(object, ...){
                         return( FittedMain(object = object@fitObj, ...) )
                       } )

setMethod(f = "ModelObjectFit", 
          signature = c(object = "QLearnEst"), 
          definition = function(object, ...){
                         return( ModelObjectFit(object = object@fitObj, ...) )
                       } )

setMethod(f = "MySummary", 
          signature = c(object = "QLearnEst"), 
          definition = function(object, ...){
                         return( MySummary(object = object@fitObj, ...) )
                       } )

setMethod(f = "Plot", 
          signature = c(x = "QLearnEst"), 
          definition = function(x, suppress=FALSE, ...){
                         Plot(x@fitObj, suppress, ...)
                       } )

setMethod(f = "PredictCont", 
          signature = c(object = "QLearnEst", 
                        newdata = "data.frame"), 
          definition = function(object, newdata, ...){
                         me <- PredictCont(object = object@fitObj, 
                                           newdata = newdata, ...)
                         return( me )
                       } )

setMethod(f = "PredictMain", 
          signature = c(object = "QLearnEst", 
                        newdata = "data.frame"), 
          definition = function(object, newdata, ...){

                         me <- PredictMain(object = object@fitObj, 
                                           newdata = newdata, ...)
                         return( me )
                       } )


setMethod(f = "Print",
          signature = c(x = "QLearnEst"),
          definition = function(x, ...){
            Print(x@fitObj)
          } )

setMethod(f = "Residuals",
          signature = c(object="QLearnEst"),
          definition = function(object, ...){ 
                         return( Residuals(object@fitObj) ) 
                       } )


setMethod(f = "Show", 
          signature = c(object = "QLearnEst"), 
          definition = function(object, ...){
                         Show(object = object@fitObj, ...)
                       } )

setMethod(f = "YTilde", 
          signature = c(object = "QLearnEst"), 
          definition = function(object, ...){ 
                         return( object@valueFunc ) 
                       } )


setMethod(f = "YTilde<-", 
          signature = c(object = "QLearnEst"), 
          definition = function(object,value){
                         object@valueFunc <- value
                         return( object )
                       } )

if(!isClass("QLearnEst or QLearnEstList")){
  setClassUnion("QLearnEst or QLearnEstList", 
                members = c("QLearnEst", "QLearnEstList")
  )
}


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                                 CLASS QLearn                                 #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# contains results obtained by fitting models for a stage of Q-learning        #
#                                                                              #
#   step       : Step in the Q-Learning algorithm                              #
#                                                                              #
#   qFunctions : matrix of Q-Functions                                         #
#                                                                              #
#   optTx      : array of optimal treatments                                   #
#                                                                              #
#   call : matched.call to QLearn function                                     #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
if(!isClass("integer or NULL")){
  setClassUnion("integer or NULL", 
                members = c("integer", "NULL")
  )
}

setClass(Class = "QLearn",
         slots = c(step = "integer",
                   qFunctions = "matrix",
                   optimalTx = "integer or factor",
                   call = "call"),
         contains = c("QLearnEst","TxInfo","DynTxRegime") )


setMethod(f = "Call",  
          signature = c(object = "QLearn"), 
          definition = function(object, ...){ 
                         return( object@call ) 
                       } )


setMethod(f = "Classif",    
          signature = c(object = "QLearn"), 
          definition = function(object, ...){
                         return( NULL )
                       } )

setMethod(f = "Est",    
          signature = c(x = "QLearn"), 
          definition = function(x, ...){
                         return( mean(x@valueFunc) )
                       } )

setMethod(f = "Genetic",    
          signature = c(object = "QLearn"), 
          definition = function(object, ...){
                         return( NULL )
                       } )

setMethod(f = "OptTx", 
          signature = c(x = "QLearn",
                        newdata = "missing"), 
          definition = function(x, newdata, ...){ 
                         return( list("qFunctions" = x@qFunctions,
                                       "optimalTx" = x@optimalTx ) )
                       } )

setMethod(f="OptTx",
          signature = c(x = "QLearn", 
                        newdata = "data.frame"),
          definition = function(x, newdata,...){
                         res <- qLearn_optTx_testSet(object = x, 
                                                     newdata = newdata)
                         return( res )
                       } )

setMethod(f = "Outcome",    
          signature = c(object = "QLearn"), 
          definition = function(object, ...){
                         return( FitObject(object = object) )
                       } )

setMethod(f = "Print",
          signature = c(x = "QLearn"),
          definition = function(x, ...){
            cat("\n")
            print(Call(x))
            cat("\n")
            Print(x@fitObj)
            cat("Mean of Value Function: ", Est(x), "\n")
          } )

setMethod(f = "Propen",    
          signature = c(object = "QLearn"), 
          definition = function(object, ...){
                         return(NULL)
                       } )

setMethod(f = "RegimeCoef",    
          signature = c(object = "QLearn"), 
          definition = function(object, ...){
                         return( NULL )
                       } )

setMethod(f = "IStep", 
          signature = c(object = "QLearn"), 
          definition = function(object, ...){ 
                         return( object@step )
                       } )

setMethod(f = "Show",
          signature = c(object = "QLearn"),
          definition = function(object, ...){
            cat("\n")
            print(Call(object))
            cat("\n")
            Print(object@fitObj)
            cat("Mean of Value Function: ", Est(object), "\n")
          } )

setMethod(f = "StdDev",    
          signature = c(object = "QLearn"), 
          definition = function(object, ...){
                         return( NULL )
                       } )


setMethod(f = "Step", 
          signature = c(object = "QLearn"), 
          definition = function(object, ...){ 
                         res <- paste("Q-Learning: step", object@step, sep=" ")
                         return( res )
                       } )


