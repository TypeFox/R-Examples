#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# Calls the residual method if specified for the solver.method                 #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setMethod(f = "residuals",
          signature = c(object="modelObjFit"),
          definition = function(object, ...){
                         tmp <- try(residuals(object@fitObj, ...), silent=TRUE)
                         if(class(tmp)=='try-error'){
                           warnMsg("residuals", class(object@fitObj))
                         } else {
                           matrix(tmp, ncol=1, dimnames=list(NULL,"residuals"))
                         }
                       })


