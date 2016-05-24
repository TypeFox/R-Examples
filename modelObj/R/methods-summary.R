#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# Calls the summary method if specified for the solver.method                  #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setMethod(f = "summary",
          signature = c(object="modelObjFit"),
          definition = function(object,...){
                         s1 <- try(summary(object@fitObj,...), silent=TRUE)
                         if(class(s1)=='try-error'){
                           warnMsg("summary", class(object@fitObj))
                         } else {
                           return(s1)
                         }
                       })


