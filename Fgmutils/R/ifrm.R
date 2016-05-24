##' @title ifrm
##' @description if the object does not exist an error will not happen.
##' @param obj the object that you want to remove
##' @param env The global environment
##' @examples
##' a = 5
##' ifrm(a)
##' ifrm(b)
##' @export
ifrm <- function(obj, env = globalenv()) {
  obj <- deparse(substitute(obj))
  if(exists(obj, envir = env)) {
    rm(list = obj, envir = env)
  }
}
