setGeneric("times",
    function(x) { standardGeneric("times")} )
setMethod(f='times',
    signature=signature(x = "x12Output"),
    definition=function(x){
      ret <- list()
      ret$original <- c(start(x@a1),end(x@a1))
      if(!all(is.na(x@forecast@estimate)))
        ret$forecast <- c(start(x@forecast@estimate),end(x@forecast@estimate))
      if(!all(is.na(x@backcast@estimate)))
        ret$backcast <- c(start(x@backcast@estimate),end(x@backcast@estimate))
      else
        ret$backcast <- NULL
      return(ret)
    })
setMethod(f='times',
    signature=signature(x = "x12Single"),
    definition=function(x)times(x@x12Output))