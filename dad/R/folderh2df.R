folderh2df <- function(foldh) {
  
  name.foldh <- as.character(match.call()$foldh)
  
  # Check of the arguments
  if (!is.folderh(foldh))
    stop(paste(name.foldh, "is not of class 'folderh'."))
  if (length(foldh) != 2)
    stop("foldh must contain two data frames.")
  
  # Change the object of class "folderh" into an object of class "folder"
  fold <- folderh2folder(foldh)
  
  # Change the object of class "folder"  into a data frame
  x <- folder2df(fold)
  colnames(x)[1] <- attr(foldh, "keys")
  
  return(x)
}
