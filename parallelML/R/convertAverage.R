#' @export
convertAverage <- function(data){
  
  # When there is only one column, apply rowMeans on that column  
  if (is.vector(data[[1]])){
  
    # Apply rowmeans on the only column
    column <- data[[1]]
    for (i in 2:length(data)){
      column <- cbind(column,data[[i]])
    }
    result <- rowMeans(column)
    
    #Return the result
    return(result)
    
    # When there is more than one column, apply rowMeans column by column
  } else {
    
    # Apply rowmeans on first column
    column <- data[[1]][,1]
    for (i in 2:length(data)){
      column <- cbind(column,data[[i]][,1])
    }
    result <- rowMeans(column)
    
    # Apply rowmeans on second column
    for (j in 2:dim(data[[1]])[2]){
      column <- data[[1]][,j]
      for (i in 2:length(data)){
        column <- cbind(column,data[[i]][,j])
      }
      result <- cbind(result,rowMeans(column))
    }
    
    # Copy the correct names
    colnames(result) <- colnames(data[[1]])
    
    # Return the result
    return(result) 
  }
}