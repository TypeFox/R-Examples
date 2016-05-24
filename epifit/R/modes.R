##' Get the type or storage mode of an object.
##'
##' This function has the almost same functionality as storage.mode except for recognition of factor. Intended for inner use in \code{showContents}, which must be able to be specified outside of \code{showContents} function.
##' @title Get the type or storage mode of an object.
##' @param x any R object.
##' @seealso \code{\link{showContents}}, \code{\link[base]{storage.mode}}.
##' @export
modes <- function(x){
  if(is.factor(x)){
    return("factor")
  } else {
    return(storage.mode(x))
  }
}
