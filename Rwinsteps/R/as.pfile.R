as.pfile <- function(x) {

  if(!is.data.frame(x) | length(dim(x)) != 2)
    stop("'x' must be a data frame with dim(x) = 2")

  varcheck(x, c("entry", "measure"))

  class(x) <- c("pfile", "data.frame")

  return(x)
}
