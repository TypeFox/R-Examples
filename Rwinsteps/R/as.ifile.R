as.ifile <- function(x) {

  if(!is.data.frame(x) | length(dim(x)) != 2)
    stop("'x' must be a data frame with dim(x) = 2")

  varcheck(x, c("entry", "measure"))

  class(x) <- c("ifile", "data.frame")

  return(x)
}
