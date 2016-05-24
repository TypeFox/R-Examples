#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                              CLASS DynTxRegime                               #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# All objects returned by methods in package DynTxRegime inherit this class.   #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setClass(Class = "DynTxRegime",
         contains = "VIRTUAL")

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                            Methods for modelObjFit                           #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# Propensity for treatments use original methods of modelObj package. These    #
# methods generalize the local methods to include objects of class modelObjFit #
# returned by package modelObj.                                                #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setMethod(f = "Coef", 
          signature = c(object = "modelObjFit"), 
          definition = function(object, ...){
                         return( coef(object = object, ...) )
                       } )

setMethod(f = "FitObject", 
          signature = c(object = "modelObjFit"), 
          definition = function(object, ...){
                         return( fitObject(object = object, ...) )
                       } )

setMethod(f = "ModelObjectFit", 
          signature = c(object = "modelObjFit"), 
          definition = function(object, ...){
                         return( object )
                       } )

setMethod(f = "MySummary", 
          signature = c(object = "modelObjFit"), 
          definition = function(object, ...){
                         return( summary(object = object, ...) )
                       } )

setMethod(f = "Plot", 
          signature = c(x = "modelObjFit"), 
          definition = function(x, suppress=FALSE, ...){
                         return( plot(x = x, ...) )
                       } )

setMethod(f = "Predict", 
          signature = c(object = "modelObjFit",
                        newdata = "data.frame"),
          definition = function(object, newdata, ...){

                         mm <- predict(object = object, newdata=newdata,...)
                         return( mm )
 
                      } )

setMethod(f = "PredictPropen", 
          signature = c(object = "modelObjFit",
                        newdata = "data.frame"),
          definition = function(object, newdata, txName, ...){

                         predArgs <- predictorArgs(object)

                         if( is(predArgs[["propen.missing"]], "NULL") ) {
                           small <- TRUE
                         } else if( predArgs[["propen.missing"]] == "smallest" ) {
                           small <- TRUE
                         } else if( predArgs[["propen.missing"]] == "largest" ) {
                           small <- FALSE
                         } else {
                           UserError("input", 
                                     "Do not recognize propen.missing value.")
                         }

                         predArgs[["propen.missing"]] <- NULL
                         predictorArgs(object) <- predArgs

                         mm <- predict(object = object, newdata=newdata,...)

                         nfv <- ncol(mm)

                         if( is(newdata[,txName], "factor") ) {
                           nlevs <- length(levels(newdata[,txName]))
                           nms <- level(newdata[,txName])
                         } else {
                           nlevs <- length(unique(newdata[,txName]))
                           nms <- as.character(sort(unique(newdata[,txName])))
                         }

                         if( nfv == {nlevs-1L} ) {

                           correction <- 1.0 - rowSums(mm)

                           if( small ) {                           
                             mm <- cbind(correction, mm)
                           } else {                           
                             mm <- cbind(mm, correction)
                           }

                         } else if( nfv != nlevs ) {
                           msg <- paste("Number of tx options in data ",
                                        "and/or fSet do not match ",
                                        "model predictions from moPropen ",
                                        "object.",sep="")
                           UserError("input", msg)
                           
                         }

                         colnames(mm) <- nms

                         return( mm )
                       } )

setMethod(f = "Print",
          signature = c(x="modelObjFit"),
          definition = function(x, ...){
                         print(x)
                       } )

setMethod(f = "Residuals", 
          signature = c(object = "modelObjFit"), 
          definition = function(object, ...){
                         return( residuals(object = object, ...) )
                       } )

setMethod(f = "Show",
          signature = c(object="modelObjFit"),
          definition = function(object, ...){
                         show(object = object)
                       } )

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                                  CLASS List                                  #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# Class mimicking an object of class list. This class is used to allow most    #
# DynTxRegime classes to be either a single object or a list of objects.       #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setClass(Class = "List",
         slots = c(loo = "list"),
         contains = "VIRTUAL")

setMethod(f = "[[",
          signature = c(x = "List"),
          definition = function(x,i){ 
                         return( x@loo[[i]] ) 
                       } )

setMethod(f = "[[<-",
          signature = c(x = "List"),
          definition = function(x,i,value){ 
                         x@loo[[i]] <- value
                         return( x ) 
                       } )

setMethod(f = "length",
          signature = c(x = "List"),
          definition = function(x){  
                         return( length(x@loo) ) 
                       } )

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                                CLASS dpList                                  #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# Class mimicking an object of class list used for multiple decision points    #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setClass(Class = "dpList",
         contains = c("List", "VIRTUAL"))

CycleThruDP <- function(object, func, ...){

  n <- length(object)

  argList <- list(...)

  res <- list()

  for( i in 1L:n ) {
    argList[[ "object" ]] <- object[[i]]
    nms <- paste("dp=", i, sep="")
    res[[ nms ]] <- do.call(what = func, args = argList)
  }

  return( res )

}

setMethod(f = "Coef", 
          signature = c(object = "dpList"), 
          definition = function(object, ...){
                         res <- CycleThruDP(object = object, 
                                            func = Coef, ...)
                         return( res )
                       } )

setMethod(f = "FitObject", 
          signature = c(object = "dpList"), 
          definition = function(object, ...){
                         res <- CycleThruDP(object = object, 
                                            func = FitObject, ...)
                         return( res )
                       } )

setMethod(f = "ModelObjectFit", 
          signature = c(object = "dpList"), 
          definition = function(object, ...){
                         res <- CycleThruDP(object = object, 
                                            func = ModelObjectFit, ...)
                         return( res )
                       } )

setMethod(f = "MySummary", 
          signature = c(object = "dpList"), 
          definition = function(object, ...){
                         res <- CycleThruDP(object = object, 
                                            func = MySummary, ...)
                         return( res )
                       } )

setMethod(f = "Plot",
          signature = c(x = "dpList"),
          definition = function(x, suppress=FALSE, ...){

                         n <- length(x)

                         for( i in 1L:n ) {

                           argList <- list(...)

                           if( !suppress ) {
                             nms <- paste("dp=", i, sep="")
                             if( is(argList[[ "main" ]], "NULL") ) {
                               argList[[ "main" ]] <- nms
                             } else if( is(argList[[ "sub" ]], "NULL") ) {
                               argList[[ "sub" ]] <- nms
                             } else {
                               argList[[ "sub" ]] <- paste(argList[[ "sub" ]], 
                                                       " (", nms, ")", sep="")
                             }
                           }
                           argList[[ "x" ]] <- x[[i]]
                           argList[[ "suppress" ]] <- suppress
                           do.call(what = Plot, args = argList)
                         }
                       } )

setMethod(f = "Residuals",
          signature = c(object = "dpList"),
          definition = function(object, ...){
                         res <- CycleThruDP(object = object, 
                                            func = Residuals, ...)
                         return( res )
                       } )

setMethod(f = "Show", 
          signature = c(object = "dpList"), 
          definition = function(object, ...){
                         n <- length(object)
                         for( i in 1L:n ) {
                           cat("Decision Point ", i, "\n")
                           Show(object = object[[i]], ...)
                         }
                       } )

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                             CLASS ModelObjList                               #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# A list of objects of class modelObj, used for multiple decision points       #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#setClass(Class = "ModelObjList",
         contains = "dpList" )



