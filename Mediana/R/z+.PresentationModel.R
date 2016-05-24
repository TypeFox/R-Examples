######################################################################################################################

# Function: +.PresentationModel.
# Argument: Two objects (PresentationModel and another object).
# Description: This function is used to add objects to the PresentationModel object
#' @export
"+.PresentationModel" = function(presentationmodel, object) {

  if (is.null(object))
    return(presentationmodel)
  else if (class(object) == "Project"){
    presentationmodel$project$username = object$username
    presentationmodel$project$title = object$title
    presentationmodel$project$description = object$description
  }
  else if (class(object) == "Section"){
    presentationmodel$section.by = unclass(object)
  }
  else if (class(object) == "Subsection"){
    presentationmodel$subsection.by = unclass(object)
  }
  else if (class(object) == "Table"){
    presentationmodel$table.by = unclass(object)
  }
  else if (class(object) == "CustomLabel"){
    ncustomlabel = length(presentationmodel$custom.label)
    presentationmodel$custom.label[[ncustomlabel+1]] = unclass(object)
  }
  else stop(paste0("Presentation Model: Impossible to add the object of class ",class(object)," to the Presentation Model"))

  return(presentationmodel)

}