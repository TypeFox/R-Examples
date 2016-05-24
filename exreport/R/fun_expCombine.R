#' Combine two experiments with different outputs
#'
#' This fuctions joints two experiments sharing the same configuration
#' of methods, problems and parameters but different outputs. The resulting
#' experiment includes the common rows for both experiments with all the
#' output columns.
#'
#' @export
#' @param e1 First experiment to combine.
#' @param e2 An second experiment to combine, must share the same 
#' config as e1.
#' @param name Optional name for the resulting experiment. If not specified
#' the new experiment will be called "e1_name U e2_name"
#' @return An new experiment with common rows and all columns.
#' 
#' 
#' @examples
#' # In this example we turn the wekaExperiment into two different experiments,
#' # with different outputs to combine them:
#' 
#' df_acc  <- wekaExperiment[,
#'            c("method", "problem", "fold", "featureSelection", "accuracy")]
#' df_time <- wekaExperiment[,
#'            c("method", "problem", "fold", "featureSelection", "trainingTime")]
#'
#' exp_acc <- expCreate(df_acc, name="acc", parameter="fold")
#' exp_time <- expCreate(df_time, name="time", parameter="fold")
#'
#' # With expCombine we can mix the two experiments:
#' expCombine(exp_acc, exp_time)
#' 
expCombine <- function(e1, e2, name=NULL){
  # PARAMETER VALIDATION:
  # Check if parameters types are correct
  if (!is.experiment(e1))
    stop(.callErrorMessage("wrongParameterError", "e1", "experiment"))
  if (!is.experiment(e2))
    stop(.callErrorMessage("wrongParameterError", "e2", "experiment"))
  if (!is.null(name) && !is.character(name))
    stop(.callErrorMessage("wrongParameterError", "name", "character or null"))
  
  # Check if the configuration for the input experiments is the same
  if(length(e1$parameters) != length(e2$parameters) 
     || !all(e1$parameters %in% e2$parameters))
    stop(.callErrorMessage("parametersDifferError"))
  
  # Any new ouptut name in e2 can have the same name as e1$method or e1$problem.
  # In that case two columns would have the same name, so the reference by name
  # to that columns would cause troubles. An error is shown in that case.
  if(any(e2$outputs %in% c(e1[[e1$method]],e1[[e1$problem]])))
    stop(.callErrorMessage("invalidOutputName","outputs of e2", "method or 
                           problem of e1"))
  
  # New parameter vector
  parameters <- e1$parameters
  
  # Merge the problems by the full configuration. We put e1$parameters in by.x
  # and by.y because this parameter takes into account the order of the elements
  # (the name of by.x will be preserved). Also, if any output have the same name
  # in e1 and e2, then the one of e2 will be renamed as that name plus ".1".
  
  commonOutputs <- e2$outputs %in% e1$outputs
  if(any(commonOutputs))
    warning(.callWarningMessage("outputsRenamed",
                                toString(e2$outputs[commonOutputs])))
  
  newDat <- merge(e1$data, e2$data, 
                  by.x=c(e1$method, e1$problem, e1$parameters), 
                  by.y=c(e2$method, e2$problem, e1$parameters),
                  suffixes = c("",".1"))
  
  # What is not a parameter, the problem or the method is an output
  namesNewDat <- colnames(newDat)
  notOutputName <- c(e1$method, e1$problem, e1$parameters)
  outputs <- namesNewDat[!(namesNewDat %in% notOutputName)]
  
  # Check if there were unmatching rows
  if (nrow(newDat) != nrow(e1$data) || nrow(newDat) != nrow(e2$data))
    warning(.callWarningMessage("unmatchedInstancesRemoved"))
  
  # Generate default name if needed
  if(is.null(name))
    newName <- paste(e1$tags$title," U ",e2$tags$title,sep="")
  else
    newName <- name
  
  # Update experiment history
  historic <- c(
    list(paste("Description for experiment ", e1$tags$title,sep="")), 
    list(e1$historic), 
    list(paste("Description for experiment ", e2$tags$title,sep="")), 
    list(e2$historic), 
    list(paste(
      "Experiment ",
      newName ,
      " created from combining two experiments: ",
      e1$tags$title," and ",e2$tags$title,".",sep="")))
  
  # Create new experiment object
  e  <- .experiment(
    data = newDat,
    method=e1$method, 
    problem=e1$problem,
    params = parameters,
    outs = outputs, 
    name = newName, 
    historic = historic)
  
  e
}
