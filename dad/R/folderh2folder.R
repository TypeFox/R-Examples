folderh2folder <- function(foldh) {        
  # Change an object foldh of class 'folderh' into an object of class 'folder',
  # Each element of this 'folder' will be the lines of x corresponding to
  # a group, and the elements of the other columns of g.
  
  name.foldh <- as.character(match.call()$foldh)
  
  # Check of the arguments
  if (!is.folderh(foldh))
    stop(paste(name.foldh, "is not of class 'folderh'."))
  if (length(foldh) != 2)
    stop("foldh must contain two data frames.")
  
  # The two data frames
  g <- foldh[[1]]
  x <- foldh[[2]]
  
  # The groups (levels of the column of g associated to key1)
  jkeyg <- which(colnames(g) == attr(foldh, "keys"))
  levg <- as.character(g[, jkeyg])
  nlevg <- length(levg)
  
  # The grouping variable in data frame x
  jkeyx <- which(colnames(x) == attr(foldh, "keys"))
  xkey <- x[, jkeyx]
  
  # Building of the list
  fold <- list()
  for (n in 1:nlevg) {
    g.tmp <- g[n, ][-jkeyg]
    fold <- c(fold, list(data.frame(x[xkey == levg[n], ][-jkeyx], g[n, ][-jkeyg], row.names = rownames(x)[xkey == levg[n]])))
  }
  names(fold) <- levg
  
  # Creation of the folder
  class(fold) <- "folder"
  attr(fold, "same.cols") <- TRUE
  attr(fold, "same.rows") <- FALSE
  
  return(fold)
}
