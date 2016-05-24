#' @export
fitted.parallelSVM <- function(object, ...){
  # Gather all the fitted results
  results <- c()
  for (i in 1:length(object)){
    column_names <- names(object[[i]]$fitted)
    temp_results <- as.character(object[[i]]$fitted)
    
    names(temp_results) <- column_names
    results      <- c(results, temp_results)
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