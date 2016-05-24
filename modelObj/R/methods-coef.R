#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
#                                                                              #
# Methods for obtaining parameter estimates from modelObjFit                   #
#                                                                              #
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::#
setMethod(f="coef",
          signature = c(object="modelObjFit"),
          definition = function(object,...){
                         ft <- fitObject(object)
                         tmp <- try(coef(ft,...), silent=TRUE)
                         if(class(tmp)=='try-error'){
                           warnMsg("coef", class(ft))
                           return(NULL)
                         } else {
                           return(tmp)
                         }
                       })




