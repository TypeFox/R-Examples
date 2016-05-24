#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
#                             CLASS ModelObjSubset                             #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# decisionPoint : decision point for which the model is to be used             #
#                                                                              #
# subset        : nickname of subset for which the model is to be used         #
#                                                                              #
# modelObject   : an object of class modelObj                                  #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#

setClass(Class = "ModelObjSubset",
         slots = c(decisionPoint = "integer",
                          subset = "character",
                     modelObject = "modelObj"))

setClass(Class = "ModelObjSubsetList",
         contains = "List" )

setMethod(f = "DecisionPoint",
          signature = c(object = "ModelObjSubset"), 
          definition = function(object, ...){ return( object@decisionPoint ) } )

setMethod(f = "Fit",
          signature = c(object = "ModelObjSubset", 
                        data = "data.frame",  
                        response = "vector"),
          definition = function(object, data, response,...){
                           res <- fit(object = object@modelObject, 
                                      data = data, 
                                      response = response)

                         return( res )

                       } )

setMethod(f = "ModelObject",
          signature = c(object = "ModelObjSubset"), 
          definition = function(object, ...){ return( object@modelObject ) } )

setMethod(f = "Subset",
          signature = c(object = "ModelObjSubset"), 
          definition = function(object, ...){ return( object@subset ) } )



