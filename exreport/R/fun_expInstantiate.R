#' Instatiate the methods in the experiment for each one of the different
#' parameter configurations.
#' 
#' When performing statistical tests or summarizing an experiment for a given
#' output variable there can be different parameter configuration for each 
#' interaction of method and problem. Once applied the desired transformations 
#' this function can be used to remove unary parameters from the experiment or to
#' instantiate the methods for each configuration.
#' 
#' If any method is instantiated the cartesian product of the method and the
#' selected parameters is performed and included in the resulting experiment as
#' the methods variable. The name of the corresponding value will indicate the
#' name of the former method and the value of each parameter instantiated.
#' 
#' @export
#' @param e The experiment object to be instantiated
#' @param parameters A vector indicating the parameters to be instantiaded.
#' If NULL or default all parameters would be considered.
#' @param removeUnary Boolean value indicating if the unary parameters will be 
#' used in an instantiation or if the column can be erased.
#' @return an experiment object
#' 
#' @examples
#' # Create an experiment from the wekaExperiment
#' experiment <- expCreate(wekaExperiment, name="test-exp", parameter="fold")
#' 
#' # We would like to reduce the fold parameter by its mean value. It becomes an
#' # unary parameter.
#' experiment <- expReduce(experiment, "fold", mean)
#' 
#' # Now we instantiate the experiment by the featureSelection parameter and
#' # remove the unary fold parameter
#' expInstantiate(experiment, removeUnary=TRUE)
#' 
expInstantiate <- function(e, parameters = NULL, removeUnary=TRUE) {
  # PARAMETER VALIDATION:
  # Check if parameters are correct
  if (!is.experiment(e))
    stop(.callErrorMessage("wrongParameterError",
                           "e",
                           "experiment"))
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
  
  # Remove unary if needed
  auxParameters <- parameters
  if (removeUnary) {
    for (var in auxParameters) {
      col <- e1$data[[var]]
      if (length(unique(col)) == 1) {
        e1$data <- e1$data[,-which(colnames(e1$data)==var)]
        parameters <- parameters[-which(parameters==var)]
      }
    }
  }
  
  # Update the data.frame
  e1$data[[e1$method]] <- interaction(e1$data[,c(e1$method,parameters),drop=FALSE],
                                      sep = ",",
                                      drop=TRUE)
  
  e1$data <- e1$data[,!(colnames(e1$data) %in% parameters),drop=FALSE]
  
  # Append this operation in the historic
  e1$historic <- c(e1$historic, 
                   list(paste("Methods has been instanciated with the parameters: ",
                              toString(auxParameters), sep="")))
 
  # Reflect this operation in the configuration to show instantiated parameters:
  config = c()
  for (p in auxParameters)
    config <- c(config, paste0(p, ' [', paste0(levels(e$data[[p]]), collapse = ","), '] (instantiated)'))
  
  e1$configuration = config
  
  e1
}

