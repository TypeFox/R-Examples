######################################################################################################################

# Function: +.DataModel.
# Argument: Two objects (DataModel and another object).
# Description: This function is used to add objects to the DataModel object
#' @export
"+.DataModel" = function(datamodel, object) {

  if (is.null(object))
    return(datamodel)

  if (class(object) == "SampleSize"){
    datamodel$general$sample.size = unclass(unlist(object))
  }
  else if (class(object) == "Event"){
    datamodel$general$event$n.events = unclass(unlist(object$n.events))
    datamodel$general$event$rando.ratio = unclass(object$rando.ratio)
  }
  else if (class(object) == "OutcomeDist"){
    datamodel$general$outcome.dist = unclass(object$outcome.dist)
    datamodel$general$outcome.type = unclass(object$outcome.type)
  }
  else if (class(object) == "Sample"){
    nsample = length(datamodel$samples)
    datamodel$samples[[nsample+1]] = unclass(object)
  }
  else if (class(object) == "Design"){
    ndesign = length(datamodel$general$design)
    datamodel$general$design[[ndesign+1]] = unclass(object)
  }
  else stop(paste0("Data Model: Impossible to add the object of class",class(object)," to the Data Model"))

  return(datamodel)

}