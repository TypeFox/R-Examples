#' Remove duplicated rows from an experiment
#'
#' This function removes duplicated rows of a given experiment attending to the
#' interaction of methods, problems and parameters (but no outputs).
#' 
#' The duplicated rows found are compared among themselves to determine if there
#' is divergence between the outputs, if the rows are not consistent a warning is
#' raised to note this difference.
#'
#' @export
#' @param e The experiment to be analised
#' @param tol The tolerance for numeric values to check if two outputs are 
#' numerically equal or not.
#' @return an experiment object
#' 
#' @examples
#' # We duplicate some of the rows of a given experiment:
#' e <- expCreate(wekaExperiment, parameters="fold", name="Test Experiment")
#' redundant <- expCreate(wekaExperiment[wekaExperiment$method=="NaiveBayes",], 
#'                        parameters="fold", name="Test Experiment")
#' e2 <- expConcat(e,redundant)
#' 
#' # Now we remove those duplicates:
#' expRemoveDuplicated(e2)
#' 
expRemoveDuplicated <- function(e, tol = 1E-9){
  # Returns a new experiment with only the duplicated instances
  #
  # Args:
  #   e:           The experiment
  #
  # Returns:
  #   a experiment object to be used in the toolbox
  
  # PARAMETER VALIDATION:
  if (!is.experiment(e))
    stop(.callErrorMessage("wrongParameterError", "e", "experiment"))
  
  # Separate data.frame in input and output part to avoid problems with numeric precission
  input <- e$data[,c(e$method,e$problem,e$parameters),drop=FALSE]
  output <- e$data[,e$outputs,drop=FALSE]
  all <- e$data
  df <- list("input" = input, "output" = output, "all" = all)
  indexDuplicates <- .duplicates(df$input)
  rowsToRemove <- c()
  if(sum(!is.na(indexDuplicates))>0){
    # There are at least one duplicated row. If there is no difference for the outputs we remove duplicated and raise a warning. 
    # If there are differences for the outputs we stop.
    # To compare outputs (because they are numeric) we use a tolerance
    # rowsToRemove is an index array, and will be used to remove not unique rows, all at once
    for(i in 1:length(indexDuplicates)){
      if(!is.na(indexDuplicates[i])){
        #The row i in df$input is the same as the row indexDuplicates[i] in df$input. Are the same outputs?
        row1 <- as.numeric(df$output[i,])
        row2 <- as.numeric(df$output[indexDuplicates[i],])
        differentOutput <- sum(abs(row1-row2)>=tol)>0
        if(differentOutput){
          stop(.callErrorMessage("failedIntegrityError"))
        }
        # If the outputs are the same, we remove this occurrence from dfs[[lesserdf]]$all
        rowsToRemove <- append(rowsToRemove, i)
      }
    }
    # Now, we remove not unique rows, all at once, and raise the warning
    df$all <- df$all[-1*rowsToRemove,,drop=FALSE]
    warning(sprintf("%d duplicated rows has been removed", length(rowsToRemove)))
  }
  result      <- e
  result$data <- df$all
  #If we removed non unique rows, we append this operation in the historic
  if(sum(!is.na(indexDuplicates))>0)
    result$historic <- c(result$historic, list(sprintf("%d duplicated rows has been removed (using a tolerance of %.4e to compare the outputs)",length(rowsToRemove),tol)))
  result
}

