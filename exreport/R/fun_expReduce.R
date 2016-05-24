#' Reduce a parameter by a function for each method, problem and remaining
#' parameter configuration interaction
#'
#' This functions reduces a parameter by aggregating the outputs variables for each
#' value and for each configuration of method, problem and remaining parameters.
#' By default it computes the mean of the variables.
#' 
#'
#' @export
#' @param e An input experiment object.
#' @param parameters The parameter or parameters to be reduced, if NULL or
#' default all parameters are considered.
#' @param FUN The function used to agregate the ouput values
#' @return An experiment object.
#' 
#' 
#' @examples
#' # Create an experiment from the wekaExperiment
#' experiment <- expCreate(wekaExperiment, name="test-exp", parameter="fold")
#' 
#' # We would like to reduce the fold parameter by its mean value. This way
#  # we compute the mean for cross validation experiment for each method and 
#  # problem, and each one of the configurations from the featureSelection param.
#' expReduce(experiment, "fold", mean)
#' 
#' 
expReduce  <- function(e, parameters = NULL, FUN=mean){
  # PARAMETER VALIDATION:
  # Check if parameters are correct
  if (!is.experiment(e))
    stop(.callErrorMessage("wrongParameterError", "e", "experiment"))
  if (!is.null(parameters) & !is.character(parameters))
    stop(.callErrorMessage("wrongParameterError",
                           "parameters",
                           "character or character array"))
  
  # If the value is NULL we set it to all the parameters
  if(is.null(parameters))
    parameters = e$parameters
  else{
    # Check if the parameters exists
    for (var in parameters) {    
      if (!(var %in% e$parameters))
        stop(.callErrorMessage("variableNotPresentError", var))
    }
  }
  
  e1 <- e
  # Update the new parameters
  e1$parameters <- e1$parameters[!(e1$parameters %in% parameters)]
  # Update the data.frame
  by <- list()
  by[[e1$method]]<-e1$data[[e1$method]]
  by[[e1$problem]]<-e1$data[[e1$problem]]
  for (p in e1$parameters) {    
    by[[p]] <- e1$data[[p]]
  }
  
  # Only apply the function to numeric values (outputs)
  extendedFUN <- function(x){
    if(is.numeric(x))
      r <- FUN(x)
    else
      r <- x
    r
  }
  # Only the aggregated columns and outputs are of interest
  e1$data <- aggregate(e1$data,
                       by=by,
                       FUN=extendedFUN)[,c(names(by),e1$outputs), drop=FALSE]
  
  #Append this operation in the historic
  e1$historic <- c(e1$historic, 
                   list(paste("Parameters '",
                              toString(parameters),
                              "' have been removed using the function '",
                              formals()$FUN,"'", sep="")))
  
  e1
}

