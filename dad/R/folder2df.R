folder2df <- function(fold) {
  
  
  name.foldh <- as.character(match.call()$foldh)
  
  # Check of the arguments
  if (!is.folder(fold))
    stop(paste(name.foldh, "is not of class 'folder'."))
  if (!attr(fold, "same.cols"))
    stop("The elements of the 'folder' must have the same number of columns and the same column names.")
  
  # The grouping variable (given by the names of the elements of the folder)
  g <- names(fold)
  
  # Building of the data frame (the grouping variable on the last column):
  # - the first group
  x <- data.frame(g = g[1], fold[[1]])
  # - and the next groups
  if (length(g) > 1) {
    for (n in 2:length(g)) {
      x <- rbind(x, data.frame(g = g[n], fold[[n]]))
    }
  }
  
  return(x)
}
