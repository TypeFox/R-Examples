df2folder <- function(x, groups = tail(colnames(x), 1)) {
  
  name.x <- as.character(match.call()$x)
  name.g <- as.character(match.call()$groups)
  
  # Checking of the arguments
  if (!is.data.frame(x))
    stop(paste(name.x, "is not a data frame."))
  if (!groups %in% colnames(x))
    stop(paste(name.g, " is not a column name of ", name.x, ".", sep = ""))
  
  # The groups
  jg <- which(colnames(x) == groups)
  g <- x[, jg]
  if (!is.factor(g)) {
    stop(paste(name.x, "[, '", name.g, "']", " is not a factor.", sep = ""))
  }
  glev <- levels(g)
  
  # Building of the list of data frames
  fold <- list()
  for (l in glev) {
    fold <- c(fold, list(x[g == l, ][-jg]))
  }
  names(fold) <- glev
  
  # Creation of the folder
  class(fold) <- "folder"
  attr(fold, "same.cols") <- TRUE
  attr(fold, "same.rows") <- FALSE
  
  return(fold)
}
