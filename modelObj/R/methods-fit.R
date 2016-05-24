#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# Obtains parameter estimates for formula objects                              #
# Returns an object of class modelObjFit                                       #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
if(!isGeneric("fit")){
  setGeneric(name = "fit", 
             def = function(object, data, response, ...){standardGeneric("fit")})
}

setMethod(f = "fit",
          signature = c(object="modelObj", data="data.frame", response="vector"),
          definition = function(object, data, response, ...){

                         #---------------------------------------------------#
                         # update formula with name of response variable     #
                         #---------------------------------------------------#
                         object@model <- update(object@model, YinternalY ~ .)

                         #---------------------------------------------------#
                         # add response to data.matrix                       #
                         #---------------------------------------------------#
                         nms <- colnames(data)
                         data <- cbind(data,response)
                         colnames(data) <- c(nms, "YinternalY")

                         #---------------------------------------------------#
                         # Set the first argument to solver as the formula   #
                         #---------------------------------------------------#
                         object@solver@methodArgs[[ 1 ]] <- object@model

                         #---------------------------------------------------#
                         # Set the data argument to the local dataset        #
                         #---------------------------------------------------#
                         object@solver@methodArgs[[ 2 ]] <- quote(data)

                         #---------------------------------------------------#
                         # Perform the fit                                   #
                         #---------------------------------------------------#
                         fit <- do.call(what=object@solver@method, 
                                        args=object@solver@methodArgs)

                         #---------------------------------------------------#
                         # Save results as new modelObjFit object and return #
                         #---------------------------------------------------#
                         ft <- new("modelObjFit",
                                   fitObj = fit, 
                                   func = object@predictor)

                         return(ft)
                       })

