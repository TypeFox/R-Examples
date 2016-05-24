#' Helper function to change default options 
#'
#' @param availableOptions Options provided by function.
#' @param optionsToSet Options to be changed to user-provided values.
adjustOptions <- function(availableOptions, optionsToSet){
  
  namesAvailableOptions <- names(availableOptions)
  changeOptions <- namesAvailableOptions[namesAvailableOptions %in% names(optionsToSet)]
  if(length(changeOptions)>0){
    for (option in changeOptions) 
      availableOptions[[option]] <- optionsToSet[[option]]
  }
  
  availableOptions
}