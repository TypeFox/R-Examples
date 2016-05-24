###############################################################################
## Functions and methods for "MCEstimate" classes
###############################################################################


setMethod("optimwarn", "MCEstimate", function(object) object@optimwarn)
setMethod("criterion", "MCEstimate", function(object) object@criterion)
setMethod("criterion.fct", "MCEstimate", function(object) object@criterion.fct)
setMethod("method", "MCEstimate", function(object) object@method)

setReplaceMethod("criterion", "MCEstimate", 
                  function(object, value) {object@criterion <- value; object})

setMethod("startPar", "MCEstimate", function(object){
 if(is.null(object@startPar)){
   cat(gettext("No starting value used.\n"))
   return(invisible(NULL))
 }else{
   return(object@startPar)
 }
 })

