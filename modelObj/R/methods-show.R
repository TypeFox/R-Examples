#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# Calls the show method if specified for the solver.method                     #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setMethod(f = "show",
          signature = c(object="modelObjFit"),
          definition = function(object){
                         show(object@fitObj)
                       })


