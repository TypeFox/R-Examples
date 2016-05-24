#' @rdname RecommenderMethods-methods
#' @aliases setMyrrixHyperParameters setMyrrixHyperParameters,list-method 
#' @exportMethod setMyrrixHyperParameters
setGeneric("setMyrrixHyperParameters", function(params, ...) standardGeneric("setMyrrixHyperParameters"))
setMethod("setMyrrixHyperParameters", signature=signature(params = "list"),
          definition = function(params){
            stopifnot(length(names(params)) == length(params) & length(list) > 0)
            javasystem <- .jnew("java.lang.System")
            for(systemproperty in names(params)){
              javasystem$setProperty(systemproperty, as.character(params[[systemproperty]]))
            }
          })

#' @rdname RecommenderMethods-methods
#' @aliases getMyrrixHyperParameters getMyrrixHyperParameters,character-method getMyrrixHyperParameters,missing-method 
#' @exportMethod getMyrrixHyperParameters 
setGeneric("getMyrrixHyperParameters", function(parameters, ...) standardGeneric("getMyrrixHyperParameters"))
.getMyrrixHyperParameters <- function(parameters){
  javasystem <- .jnew("java.lang.System")            
  systemproperties <- list()
  props <- javasystem$getProperties()$entrySet()$iterator()
  while(props$hasNext()){
    entry <- props$nextElement()
    systemproperties[[entry$getKey()]] <- entry$getValue()
  }           
  systemproperties
}
setMethod("getMyrrixHyperParameters", signature="missing", definition = .getMyrrixHyperParameters)
setMethod("getMyrrixHyperParameters", signature="character", definition = function(parameters){
  systemproperties <- .getMyrrixHyperParameters(parameters)
  if(length(parameters) > 0){
    systemproperties <- systemproperties[names(systemproperties)[names(systemproperties) %in% parameters]]
  } 
  systemproperties
})

