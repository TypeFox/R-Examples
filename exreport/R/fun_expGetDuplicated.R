#' Create a new experiment with only the duplicated rows
#'
#' This function computes the duplicated rows attending to the method, problem
#' and input parameters (but not the outputs). The resulting experiment will
#' contain these duplicated rows. 
#' 
#' If duplicated rows show different outputs the function will launch a 
#' a warning message indicating how many of them differ in the outputs from the
#' original row, the extent to what two rows are divergent in their output can be
#' parametrized. 
#' 
#' This function is useful to determine the consistency of the experiment, as a
#' measure to sanitice the original data source if needed,
#'
#' @export
#' @param e The experiment to check for duplicated rows
#' @param tol The tolerance for numeric values to check if two outputs are 
#' numerically equal or not.
#' @return A new experiment containing the duplicated rows
#' 
#' 
#' @examples
#' # We duplicate some of the rows of a given experiment:
#' e <- expCreate(wekaExperiment, parameters="fold", name="Test Experiment")
#' redundant <- expCreate(wekaExperiment[wekaExperiment$method=="NaiveBayes",], 
#'                        parameters="fold", name="Test Experiment")
#' e2 <- expConcat(e,redundant)
#' 
#' # Now we check for duplicates:
#' expGetDuplicated(e2)
#' 
expGetDuplicated <- function(e, tol = 1E-9){
  # PARAMETER VALIDATION:
  if (!is.experiment(e))
    stop(.callErrorMessage("wrongParameterError", "e", "experiment"))
  
  # Separate data.frame in input and output part to avoid problems with 
  # numeric precission
  input  <- e$data[,c(e$method,e$problem,e$parameters),drop=FALSE]
  output <- e$data[,e$outputs,drop=FALSE]
  all    <- e$data
  df <- list("input" = input, "output" = output, "all" = all)
  
  # The function .duplicates returns the index of rows which are duplicated
  indexDuplicates <- .duplicates(df$input)
  duplicatedRows <- c()
  diffOutputs <- 0
  if(sum(!is.na(indexDuplicates))>0){
    # There are at least one duplicated row. If there is no difference for the 
    # outputs we remove duplicated and raise a warning. 
    # If there are differences for the outputs we stop.
    # To compare outputs (because they are numeric) we use a tolerance
    # duplicatedRows is an index array, and will be used to remove not 
    # unique rows, all at once
    for(i in 1:length(indexDuplicates)){
      if(!is.na(indexDuplicates[i])){
        # The row i in df$input is the same as the row indexDuplicates[i]
        # in df$input. Are the same outputs?
        row1 <- as.numeric(df$output[i,])
        row2 <- as.numeric(df$output[indexDuplicates[i],])
        differentOutput <- sum(abs(row1-row2)>=tol)>0
        if(differentOutput){
          diffOutputs <- diffOutputs+1
        }
        # If the outputs are the same, we remove this occurrence 
        # from dfs[[lesserdf]]$all
        duplicatedRows <- append(duplicatedRows, i)
      }
    }
  }
  # Now, we create the data.frame with the duplicates
  originalOfDuplicated <- unique(indexDuplicates[!is.na(indexDuplicates)])
  df$all <- df$all[c(originalOfDuplicated,duplicatedRows),,drop=FALSE]
  
  result      <- e
  result$data <- df$all
  #If we removed non unique rows, we append this operation in the historic
  if(length(duplicatedRows)>0){
    message = sprintf("%d duplicated rows. 
                      %d differ in the outputs (using a tolerance of %.4e to compare the outputs)",length(duplicatedRows), diffOutputs,tol)
    if(diffOutputs>0)
      warning(paste0(message,". Please check the integrity of the experiment"))
    else
      warning(message)
    result$historic <- c(result$historic, list(paste0("Subset with only duplicated rows. ",message)))
  }
  result
}