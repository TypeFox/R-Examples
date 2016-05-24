#' Extend an experiment by adding new parameters
#'
#' This function extends an existing exreport experiment object by adding new 
#' parameters with fixed values.
#'
#' @export
#' @param e Input experiment
#' @param parameters A list of strings containing the values of the new 
#' parameters, the name for each one of them will be given by the name of the 
#' corresponding object in the list.
#' @return A modified exreport experiment object with additional parameters.
#'
#' @examples
#' # We load the wekaExperiment problem as an experiment and then add a new param
#' # with a default value.
#' 
#' experiment <- expCreate(wekaExperiment, name="test", parameter="fold")
#' expExtend(experiment, list(discretization = "no"))
#' 
expExtend  <- function(e, parameters){
  # PARAMETER VALIDATION:
  # Check if parameters are correct
  if (!is.experiment(e))
    stop(.callErrorMessage("wrongParameterError", "e", "experiment"))
  if (!is.list(parameters))
    stop(.callErrorMessage("wrongParameterError", "parameters", "non-empty list"))
  if (length(parameters)==0)
    stop(.callErrorMessage("wrongParameterError", "parameters", "non-empty list"))
  # Check that all parameters have a proper name
  if (is.null(names(parameters)))
    stop(.callErrorMessage("noNamesError"))
  
  # Copy the experiment
  result <- e
  # Modify the parameter vector of the experiment
  pnames <- names(parameters)
  result$parameters <- c(result$parameters, pnames)
  nrows <- nrow(result$data)
  
  # Add the new columns to the data.frame
  for (i in 1:length(parameters)){
    name <- pnames[i]
    value <- parameters[[i]]
    result$data[[name]] <- as.factor(rep(value,nrows))
  }
  
  # Append this operation in the historic
  pairNames <- paste(pnames,parameters,sep=":")
  result$historic <- c(result$historic, 
                       list(paste("New parameters have been added with default values: ",
                                  toString(pairNames) ,sep="")))
  
  result
}