#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                               CLASS SubsetFit                                #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setClass(Class = "SubsetFit",
         slots = c(  subset = "character",
                     fitObj = "SimpleFit or IterateFit"))

setClass(Class = "SubsetFitList",
         slots = c(decisionPoint = "integer",
                          txInfo = "TxInfo"),
         contains = "List" )

setMethod(f = "Coef",
          signature = c(object = "SubsetFit"),
          definition = function(object, ...) { 
                         return( Coef(object = object@fitObj, ...) ) 
                       } )

setMethod(f = "FitObject", 
          signature = c(object = "SubsetFit"), 
          definition = function(object, ...){
                         return( FitObject(object = object@fitObj, ...) )
                       } )

setMethod(f = "FittedCont", 
          signature = c(object = "SubsetFit"), 
          definition = function(object, ...){
                           return( FittedCont(object = object@fitObj, ...) )
                       } )

setMethod(f = "FittedMain", 
          signature = c(object = "SubsetFit"), 
          definition = function(object, ...){
                         return( FittedMain(object = object@fitObj, ...) )
                       } )

setMethod(f = "ModelObjectFit", 
          signature = c(object = "SubsetFit"), 
          definition = function(object, ...){
                         return( ModelObjectFit(object = object@fitObj, ...) )
                       } )

setMethod(f = "MySummary",
          signature = c(object = "SubsetFit"),
          definition = function(object, ...) {
                         return( MySummary(object = object@fitObj, ...) )
                       } )

setMethod(f = "Plot", 
          signature = c(x = "SubsetFit"), 
          definition = function(x, suppress=FALSE, ...){
                         Plot(x@fitObj, suppress, ...)
                       } )

setMethod(f = "PredictCont", 
          signature = c(object = "SubsetFit", 
                        newdata = "data.frame"), 
          definition = function(object, newdata, ...){

                         res <- PredictCont(object = object@fitObj, 
                                            newdata = newdata, ...)
                         return( res )
                       } )

setMethod(f = "PredictMain", 
          signature = c(object = "SubsetFit", 
                        newdata = "data.frame"), 
          definition = function(object, newdata, ...){

                         res <- PredictMain(object = object@fitObj, 
                                            newdata = newdata, ...)
                         return( res )
                       } )

setMethod(f = "Print", 
          signature = c(x = "SubsetFit"), 
          definition = function(x, suppress=FALSE, ...){
                         cat("Subset: ", Subset(x), "\n")
                         Print(x@fitObj)
                       } )

setMethod(f = "Residuals",
          signature = c(object = "SubsetFit"),
          definition = function(object, ...) {
                         return( Residuals(object = object@fitObj) )
                       } )

setMethod(f = "Show", 
          signature = c(object = "SubsetFit"), 
          definition = function(object, ...){
                         cat("Subset: ", Subset(object), "\n")
                         Show(object = object@fitObj, ...)
                       } )

setMethod(f = "Subset",
          signature = c(object = "SubsetFit"), 
          definition = function(object, ...){ 
                         return( object@subset ) 
                       } )


CycleThruTx <- function(object, func, ...){

  n <- length(object)

  argList <- list(...)

  res <- list()

  for( i in 1L:n ) {
    argList[[ "object" ]] <- object[[i]]
    nms <- paste(Subset(object[[i]]), collapse=",")
    res[[ nms ]] <- do.call(what = func, args = argList)
  }

  return( res )
}


setMethod(f = "Coef",
          signature = c(object = "SubsetFitList"),
          definition = function(object, ...){
                         res <- CycleThruTx(object = object, func = Coef, ...)
                         return( res )
                       } )

setMethod(f = "DecisionPoint",
          signature = c(object = "SubsetFitList"), 
          definition = function(object, ...){ 
                         return( object@decisionPoint ) 
                       } )

setMethod(f = "FitObject",
          signature = c(object = "SubsetFitList"),
          definition = function(object, ...){
                         res <- CycleThruTx(object = object, func = FitObject, ...)
                         return( res )
                       } )

setMethod(f = "FittedCont", 
          signature = c(object = "SubsetFitList"), 
          definition = function(object, ...){

                         lst <- CycleThruTx(object = object, 
                                            func = FittedCont, ...)
                         n <- length(object)

                         ptsSubset <- PtsSubset(object@txInfo)
                         res <- numeric(length(ptsSubset))

                         for( i in 1L:n ) {
                           u4f <- ptsSubset %in% Subset(object[[i]])
                           res[u4f] <- lst[[i]]
                         }

                         return( res )
                       } )

setMethod(f = "FittedMain", 
          signature = c(object = "SubsetFitList"), 
          definition = function(object, ...){

                         lst <- CycleThruTx(object = object, 
                                            func = FittedMain, ...)
                         n <- length(object)
                         ptsSubset <- PtsSubset(object@txInfo)
                         res <- numeric(length(ptsSubset))

                         for( i in 1L:n ) {
                           u4f <- ptsSubset %in% Subset(object[[i]])
                           res[u4f] <- lst[[i]]
                         }

                         return( res )

                       } )

LocalPredictFunction <- function(object, 
                                 newdata, 
                                 func, 
                                 ...){

  fs <- feasibility(superSet = SuperSet(object@txInfo), 
                    fSet = SubsetRule(object@txInfo), 
                    txName = TxName(object@txInfo),
                    data = newdata)

  ptSubsets <- fs$subsets

  nPtSubsets <- length(ptSubsets)

  res <- numeric(nrow(newdata))
  res[] <- NA

  for( i in 1L:nPtSubsets ) {

    namePtSubset <- names(ptSubsets)[i]

    iss <- 0L
    for( j in 1L:length(object) ){
      fitSubset <- Subset(object[[j]])
      if( namePtSubset %in% fitSubset ) {
        iss <- j
        break
      }
    }

    if(iss == 0L) stop("Could not match subset identified in new data.")

    u4f <- fs$ptsSubset == namePtSubset 

    if( sum(u4f) == 0L ) next

    argList <- list(...)
    argList[[ "object" ]] <- object[[j]]
    argList[[ "newdata" ]] <- newdata[u4f,,drop=FALSE]

    cn <- do.call(what = func, args = argList)

    res[u4f] <- cn
  }
  return(res)
}

LocalResidualsFunction <- function(object, ...){

  ptsSubsets <- PtsSubset(object@txInfo)

  res <- numeric(length(ptsSubsets))
  res[] <- NA

  for( i in 1L:length(object) ) {

    u4f <- ptsSubsets %in% Subset(object[[i]])

    if( sum(u4f) == 0L ) next

    cn <- Residuals(object[[i]])

    res[u4f] <- cn
  }
  return(res)
}

setMethod(f = "MySummary",
          signature = c(object = "SubsetFitList"),
          definition = function(object, ...){
                         res <- CycleThruTx(object = object, 
                                            func = MySummary, ...)
                         return( res )
                       } )

setMethod(f = "Plot",
          signature = c(x = "SubsetFitList"),
          definition = function(x, suppress=FALSE, ...){
                         n <- length(x)
                         for( i in 1L:n ) {
                           argList <- list(...)
                           if( !suppress ) {
                             nms <- paste("subset={",
                                          paste(Subset(x[[i]]),collapse=","),
                                          "}",sep="")
                             if( is(argList[[ "main" ]], "NULL") ) {
                               argList[[ "main" ]] <- nms
                             } else if( is(argList[[ "sub" ]], "NULL") ) {
                               argList[[ "sub" ]] <- nms
                             } else {
                               argList[[ "sub" ]] <- paste(argList[[ "sub" ]], 
                                                     " (",nms,")", sep="")
                             }
                           }
                           argList[[ "x" ]] <- x[[i]]
                           argList[[ "suppress" ]] <- suppress
                           do.call(what = Plot, args = argList)
                         }
                       } )

setMethod(f = "PredictCont", 
          signature = c(object = "SubsetFitList", 
                        newdata = "data.frame"), 
          definition = function(object, newdata, ...){

                         me <- LocalPredictFunction(object = object, 
                                                    newdata = newdata, 
                                                    func = PredictCont, ...)

                         return(me)
                       } )

setMethod(f = "PredictMain", 
          signature = c(object = "SubsetFitList", 
                        newdata = "data.frame"), 
          definition = function(object, newdata, ...){

                         me <- LocalPredictFunction(object = object, 
                                                    newdata = newdata, 
                                                    func = PredictMain, ...)

                         return(me)
                       } )

setMethod(f = "Print",
          signature = c(x = "SubsetFitList"),
          definition = function(x, ...){

                         n <- length(x)

                         for( i in 1L:n ) {
                           Print(x[[i]])
                         }
                         return
                       } )


setMethod(f = "Residuals",
          signature = c(object = "SubsetFitList"),
          definition = function(object, ...){
                         res <- LocalResidualsFunction(object =  object, ...)
                         return( res )
                       } )

setMethod(f = "Show",
          signature = c(object = "SubsetFitList"),
          definition = function(object, ...){
                         n <- length(object)
                         for( i in 1L:n ) {
                           Show(object[[i]], ...)
                         }
                       } )

if(!isClass("SimpleFit, IterateFit, SubsetFitList")){
  setClassUnion("SimpleFit, IterateFit, SubsetFitList", 
                members = c("SimpleFit", "IterateFit", "SubsetFitList"))
}


