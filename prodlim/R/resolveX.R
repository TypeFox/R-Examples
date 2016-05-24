resolveX <- function(object,N){
  if (missing(object)) X <- NULL
  if (!missing(object) && (is.null(object)|| (is.logical(object) && object==FALSE)))
    X <- NULL
  else{
    ## if the object is a matrix then do nothing
    if (is.matrix(object) && NROW(object)==N)
      X <- object
    else
      X <- data.frame(sapply(object, function(x) {
        ## each entry is either a distribution to draw from
        if (is.character(x[[1]]) || is.function(x[[1]]))
          do.call(x[[1]], c(n = N, x[-1]))
        else{
          ## or a vector of numeric values
          stopifnot(is.numeric(x) && length(x)==N)
          x}
      }))
  }
  X
}
