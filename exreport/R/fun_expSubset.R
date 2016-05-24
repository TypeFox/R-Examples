#' Obtains a subset of an experiment matching the given conditions
#'
#' This function receives a named list indicating variables and values to filter
#' the input experiment.
#' 
#' The names of the elements in the list correspond with the variables to be
#' filtered, indicating either the methos or problem variables as well as 
#' parameters. The values of the list correspond with the valid states for the
#' filtering.
#'
#' @export
#' @param e The experiment to be subsetted
#' @param columns A named list containing the variables to be filtered and the
#' valid values.
#' @param invertSelection If the filtering must match the inversion of the 
#' specified conditions.
#' @return a filtered experiment object
#' 
#' @examples
#' # We create a new experiment from the wekaExperiment problem
#' e <- expCreate(wekaExperiment, parameters="fold", name="Test Experiment")
#' 
#' # We can filter the experiment to reduce the number of methods.
#' e <- expSubset(e, list(method = c("J48", "NaiveBayes")))
#' e
#' 
#' # We can filter the experiment to remove a given problem
#' e <- expSubset(e, list(problem = "iris"), invertSelection=TRUE)
#' e
#' 
#' # We can subset the experiment to obtain a specific parameter configuration
#' e <- expSubset(e, list("featureSelection" = "no"))
#' e
#' 
expSubset <- function(e, columns, invertSelection=FALSE){

  
  # PARAMETER VALIDATION:
  # Check if parameters are correct
  if (!is.experiment(e))
    stop(.callErrorMessage("wrongParameterError", "e", "experiment"))
  if (!is.list(columns) | length(columns)==0)
    stop(.callErrorMessage("wrongParameterError", "columns", "non-empty list"))
  if (is.null(names(columns)))
    stop(.callErrorMessage("noNamesError"))
  
  # Check if the columns exists
  for (var in names(columns)) {    
    if (var!=e$method && var!=e$problem && !(var %in% e$parameters))
      stop(.callErrorMessage("variableNotPresentError", var))
  }
  
  result <- e
  pnames <- names(columns)
  nrows <- nrow(result$data)
  valueNames <- c()
  # First we generate an array of logical values to make the subsetting
  subset <- rep(TRUE,nrows)
  for (i in 1:length(columns)){
    name <- pnames[i]
    value <- columns[[i]]
    valueNames <- append(valueNames,paste("[",toString(value),"]", sep=""))
    subset <- subset & result$data[[name]] %in% value
  }
  
  if (invertSelection)
    subset <- !subset
  
  # Now we subset the new columns to the data.frame
  result$data <- result$data[subset,,drop=FALSE]
  
  # We append this operation in the historic
  pairNames <- paste(pnames,valueNames,sep=":")
  result$historic <- c(result$historic, list(paste("A set of pairs 'column:[defaultValues]' has been used to subset the experiment: ", toString(pairNames) ,sep="")))
  
  result
}