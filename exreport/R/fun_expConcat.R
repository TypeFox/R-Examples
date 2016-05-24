#' Concatenate rows of matching experiments
#'
#' This function concatenates two experiments with the same configuration
#' of parameter an outputs. At least one common output must be present,
#' the rest of them will be removed from the resulting experiment.
#' Different methods and problems can be present.
#'
#' @export
#' @param e1 First experiment object to concat.
#' @param e2 Second experiment object to concat. Must have the same configuration
#' than e1.
#' @param name Optional name, if not provided the new experiment will be called 
#' "e1_name + e2_name"
#' @param tol Tolerance value for duplicate checking.
#' @return An experiment object having all the rows of e1 and e2
#' 
#' 
#' @examples
#' # In this example we turn the wekaExperiment into two different experiments,
#' # with different parameter values to combine them:
#' 
#' df_no  <- wekaExperiment[wekaExperiment$featureSelection=="no",]
#' df_yes <- wekaExperiment[wekaExperiment$featureSelection=="yes",]
#'
#' exp_yes <- expCreate(df_yes, name="fss-yes", parameter="fold")
#' exp_no <- expCreate(df_no, name="fss-no", parameter="fold")
#' 
#' expConcat(exp_yes, exp_no)
#' 
expConcat <- function(e1, e2, name=NULL, tol = 1E-9){
  # PARAMETER VALIDATION:
  # Check if parameters are correct
  if (!is.experiment(e1))
    stop(.callErrorMessage("wrongParameterError", "e1", "experiment"))
  if (!is.experiment(e2))
    stop(.callErrorMessage("wrongParameterError", "e2", "experiment"))
  if (!is.null(name) && !is.character(name))
    stop(.callErrorMessage("wrongParameterError", "name", "character or null"))
  
  # Check if configurations match
  if( length(e1$parameters) != length(e2$parameters) 
      || !all(e1$parameters %in% e2$parameters))
    stop(.callErrorMessage("parametersDifferError"))
  
  # Check at least one common output is present in both experiments
  intersection <- intersect(e1$outputs, e2$outputs)
  if(length(intersection)==0)
    stop(.callErrorMessage("commonOutputsError"))
  
  # Remove uncommon outputs (if it is the case we raise a warning).
  data1           <- e1$data
  subsetE1        <- data1[,c(e1$method,e1$problem,e1$parameters,intersection), drop=FALSE]
  data2           <- e2$data
  subsetE2        <- data2[,c(e2$method,e2$problem,e2$parameters,intersection), drop=FALSE]
  if( ncol(subsetE1)!= ncol(e1$data) ||  ncol(subsetE2)!= ncol(e2$data)) {
    data1 <- subsetE1
    data2 <- subsetE2
    warning(.callWarningMessage("outputsRemoved",toString(intersection)))
  }
  
  # If the method and/or problem column names are different for e1 and e2, we stop!
  if(e1$method != e2$method){
    colnames(data2)[colnames(data2)==e2$method] <- e1$method
    warning(.callWarningMessage("differentMethodName","e1","e2",e1$method))
  }
  if(e1$problem != e2$problem){
    colnames(data2)[colnames(data2)==e2$problem] <- e1$problem
    warning(.callWarningMessage("differentProblemName","e1","e2",e1$problem))
  }
  
  # Append both data.frames
  appendedDataFrame <- rbind(data1,data2)  
  parameters <- e1$parameters
  
  # Generate new name if needed
  if(is.null(name))
    newName <- paste(e1$tags$title,"+",e2$tags$title,sep="")
  else
    newName <- name
  
  # Update experiment history
  historic <- c(list(paste("Description for experiment ", e1$tags$title,sep="")),
                list(e1$historic), 
                list(paste("Description for experiment ", e2$tags$title,sep="")), 
                list(e2$historic), 
                list(
                  paste("Experiment ",
                        newName,
                        " created from appending two experiments: ",
                        e1$tags$title," and ",e2$tags$title,".",sep="")))
  
  e  <- .experiment(
    data = appendedDataFrame,
    method=e1$method,
    problem=e1$problem,
    params = parameters,
    outs = intersection, 
    name = newName,
    historic = historic)
  
  # Check if rows are duplicated, but do not revome them.
  expGetDuplicated(e,tol=tol)
  
  e
}