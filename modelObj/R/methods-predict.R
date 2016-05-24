#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# Predict - returns matrix of predicted values obtained from a modelObjFit     #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

setMethod(f = "predict",
          signature = c(object="modelObjFit"),
          definition = function(object, newdata=NULL, ...){
                         if(!is.null(newdata) && !is.data.frame(newdata)) {
                           stop("newdata must be a data.frame")
                         }
                         object@func@methodArgs[[ 1 ]] <- object@fitObj
                         object@func@methodArgs[[ 2 ]] <- newdata
                         mm <- do.call(what = object@func@method,
                                       args = object@func@methodArgs)

                         if(!is.matrix(mm)) mm <- matrix(mm,ncol=1)

                         return(mm)
                       })


