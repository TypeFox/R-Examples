#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# valueFuncs - Calculates/Retrieves the value functions for each treatment     #
#              option for IQ- and Q-Learning methods.                          #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
if(!isGeneric("qFuncs")){
  setGeneric(name = "qFuncs", 
             def = function(object, ...){standardGeneric("qFuncs")})
}

setMethod(f = "qFuncs",
          signature = c(object="IQEst"),
          definition = function(object, ...){

                         opts <- object@qFunctions
                         colnames(opts) <- c("-1","1")

                         return(opts)
                       } )

setMethod(f = "qFuncs",
          signature = c(object="IQLearnFS_VHom"),
          definition = function(object, ...){
                         return(StdDev(object))
                        } )

setMethod(f = "qFuncs",
          signature = c(object="QLearnEst"),
          definition = function(object, ...){

                         return(object@qFunctions)
                       } )


setMethod(f = "qFuncs",
          signature = c(object="OptimalSeq"),
          definition = function(object, ...){
                         return(NULL)
                       } )

setMethod(f = "qFuncs",
          signature = c(object="OptimalClass"),
          definition = function(object, ...){
                         return(NULL)
                       } )


