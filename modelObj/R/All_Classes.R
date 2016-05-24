if(!isClass("character or NULL")){
  setClassUnion("character or NULL", members = c("character","NULL"))
}
if(!isClass("numeric or NULL")){
  setClassUnion("numeric or NULL", members = c("numeric","NULL"))
}

################################################################################
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Class methodObj and its methods                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# methodObj contains a function name and a list of arguments to be passed to   #
# the function on call.                                                        #
#                                                                              #
#   method - function name or character string containing the function name    #
#                                                                              #
#   methodArgs - list of formals of method to be passed to method when called  #
################################################################################
setClass("methodObj", 
         slots = c(    method = "character",
                   methodArgs = "list") )

setGeneric(name = "method",
           def = function(object,...){standardGeneric("method")})
setGeneric(name = "methodArgs",
           def = function(object,...){standardGeneric("methodArgs")})

setMethod(f = "method", 
          signature = c(object="methodObj"), 
          definition = function(object,...){object@method})
setMethod(f = "methodArgs", 
          signature = c(object="methodObj"), 
          definition = function(object,...){object@methodArgs})

setValidity(Class = "methodObj", 
            method = function(object){
                       if( method(object) == 'predict' ) return(TRUE)
                       nms <- names( methodArgs(object) )
                       fmls <- names( formals( method(object) ) )
                       if( !exists(method(object)) ) {
                         return("Method failed exists() check.")
                       }
                       if( all(nms %in% fmls) ) {
                         return(TRUE)
                       } else {
                         print(nms)
                         print(fmls)
                         msg <- paste("Elements of methodArgs not found in ",
                                      "formal arguments of method",sep="")
                         return(msg)
                       }
                     })
                       

################################################################################
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# modelObjFit class                                                            #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# fit is a list whose first element contains:                                  #
#  contents depend on the solverMethod specified.                              #
#                                                                              #
# func is the predictMethod corresponding to solverMethod                      #
################################################################################
setClass("modelObjFit",
         slots = c(fitObj = "ANY",
                     func = "methodObj"))

setGeneric(name = "fitObject",
           def = function(object,...){standardGeneric("fitObject")})
setGeneric(name = "predictor",
           def = function(object,...){standardGeneric("predictor")})
setGeneric(name = "predictorArgs",
           def = function(object,...){standardGeneric("predictorArgs")})

setMethod(f = "fitObject", 
          signature = c(object="modelObjFit"), 
          definition = function(object,...){object@fitObj})
setMethod(f = "predictor", 
          signature = c(object="modelObjFit"), 
          definition = function(object,...){object@func@method})
setMethod(f = "predictorArgs", 
          signature = c(object="modelObjFit"), 
          definition = function(object,...){object@func@methodArgs})

################################################################################
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# Class modelObj                                                               #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
# contains methods and arguments for fitting                                   #
#                                                                              #
#  model : a single formula objects to be evaluated                            #
#                                                                              #
#  solver : methodObj for parameter estimates                                  #
#                                                                              #
#  predictor : methodObj for predictions                                       #
################################################################################
setClass("modelObj",
    slots = c(    model = "formula",
                      solver = "methodObj",
                   predictor = "methodObj")
)
setGeneric(name = "model",
           def = function(object,...){standardGeneric("model")})
setGeneric(name = "solver",
           def = function(object,...){standardGeneric("solver")})
setGeneric(name = "solverArgs",
           def = function(object,...){standardGeneric("solverArgs")})

setGeneric(name = "solverArgs<-",
           def = function(object,value){standardGeneric("solverArgs<-")})
setGeneric(name = "predictorArgs<-",
           def = function(object, value){standardGeneric("predictorArgs<-")})

setMethod(f = "model", 
          signature = c(object="modelObj"), 
          definition = function(object,...){object@model})
setMethod(f = "solver",  
          signature = c(object="modelObj"), 
          definition = function(object,...){object@solver@method})
setMethod(f = "solverArgs",  
          signature = c(object="modelObj"), 
          definition = function(object,...){object@solver@methodArgs})
setMethod(f = "predictor",  
          signature = c(object="modelObj"), 
          definition = function(object,...){object@predictor@method})
setMethod(f = "predictorArgs",  
          signature = c(object="modelObj"), 
          definition = function(object,...){object@predictor@methodArgs})

setMethod(f = "solverArgs<-",   
          signature = c(object="modelObj"), 
          definition = function(object, value){
                         nms <- names(object@solver@methodArgs)
                         nmsNew <- names(value)
                         if(!(nms[2] %in% nmsNew)) {
                           value <- c(object@solver@methodArgs[[2]],value)
                         }
                         if(!(nms[1] %in% nmsNew)) {
                           value <- c(object@solver@methodArgs[[1]],value)
                         }
                         object@solver@methodArgs <- value
                         return(object)
                       })


setMethod(f = "predictorArgs<-",   
          signature = c(object="modelObj"), 
          definition = function(object, value){
                         nms <- names(object@predictor@methodArgs)
                         nmsNew <- names(value)
                         if(!(nms[2] %in% nmsNew)) {
                           value <- c(object@predictor@methodArgs[[2]],value)
                         }
                         if(!(nms[1] %in% nmsNew)) {
                           value <- c(object@predictor@methodArgs[[1]],value)
                         }
                         object@predictor@methodArgs <- value
                         return(object)
                       })

