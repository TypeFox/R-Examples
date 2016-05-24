#' bmr class
#' 
#' Creates an object of class bmr (bootstrap metaregression)
#' 
#'  @section Slots:
#'    \describe{
#'      \item{\code{models}:}{Object of class \code{"list"}, containing the
#'      models.}
#'      \item{\code{fits}:}{Object of class \code{"list"}, containing the model
#'      fit predictions.}
#'      \item{\code{medians}:}{Object of class \code{"vector"}, containing the
#'      medians.}
#'      \item{\code{confidence_envelope}:}{Object of class \code{"data.frame"}, 
#'      containing the lower and upper confidence bounds.}
#'    }
#'  @name bmr-class
#'  @rdname bmr-class
#'  @aliases bmr
#'  @exportClass bmr
#'  @author Lyle D. Burgoon
#'  
#'  
#'  @export
bmr <- setClass(
  #Set the name
  "bmr",
  
  #Define the slots
  slots = c(
    models = "list",
    fits = "list",
    medians = "vector",
    confidence_envelope = "data.frame"
  ),
  
  validity=function(object){
    if(!is.list(object@models)){
      return("models must be a list")
    }
    if(!is.list(object@fits)){
      return("fits must be a list")
    }
    if(!is.vector(object@medians)){
      return("medians must be a vector")
    }
    if(!is.data.frame(object@confidence_envelope)){
      return("confidence envelope must be a data frame")
    }
    return(TRUE)
  }
  
)