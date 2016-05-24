#' @export
fitted.parallelML <- function(object, ...){
  # Gather all the fitted results
  results <- c()
  
  if ("fitted" %in% names(object[[1]])){
    for (i in 1:length(object)){
      column_names <- names(object[[i]]$fitted)
      temp_results <- as.character(object[[i]]$fitted)
      
      names(temp_results) <- column_names
      results             <- c(results, temp_results)
    }
    warning("Not all observations might be fitted, due to the sampling.")
  } else if ("fitted.values" %in% names(object[[1]])){
    if (!is.null(names(object[[1]]$fitted.values))){
      for (i in 1:length(object)){
        column_names <- names(object[[i]]$fitted.values)
        temp_results <- as.character(object[[i]]$fitted.values)
        
        names(temp_results) <- column_names
        results             <- c(results, temp_results)
      }
      warning("Not all observations might be fitted, due to the sampling.")
    } else {
      for (i in 1:length(object)){
        column_names <- as.numeric(attr(object[[i]]$fitted.values,"dimnames")[[1]])
        temp_results <- as.character(object[[i]]$fitted.values)
        
        names(temp_results) <- column_names
        results             <- c(results, temp_results)
      }
    }
  } else {
    stop("No fitted values found in your model")
  }
  
  # Delete .1 from names
  column_names         <- round(as.numeric(names(results)))
  
  # Convert to a useful table
  results              <- data.frame(results)
  results$column_names <- column_names
  results              <- table(results)
  
  # Get most occuring lable
  final_result <- c()
  column_names <- colnames(results)
  for (i in 1:length(column_names)){
    name <- column_names[i]
    final_result[i] <- names(sort(results[,name],decreasing=TRUE))[1]
  }
  names(final_result) <- column_names
  
  final_result <- as.factor(final_result)
  
  return(final_result)
}