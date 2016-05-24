#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# Calls the show method if specified for the solver.method                     #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setMethod(f = "print",
          signature = c(x="modelObjFit"),
          definition = function(x){
                         print(x@fitObj)
                       })


