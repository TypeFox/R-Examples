listn <- function(...) {
  input <- sapply(match.call(), deparse)[-1]
  if (is.null(names(input))) {
    names(input) <- rep("", times=length(input))
  }
  names(input) <- ifelse(names(input)=="", input, names(input))
  listWithName <- list(...)
  names(listWithName) <- names(input) 
  return(listWithName)
}