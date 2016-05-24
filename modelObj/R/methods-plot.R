#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# Calls the plot method if specified for the solver.method                     #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setMethod(f = "plot",
          signature = c(x="modelObjFit"),
          definition = function(x, ...){
                         tmp <- try( plot(x@fitObj, ...), silent=TRUE )
                         if( class(tmp) == 'try-error' ){
                           warnMsg("plot", class(x@fitObj))
                         }
                       })


