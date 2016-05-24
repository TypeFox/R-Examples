chklist <- function(x){
  if (!is.list(x))
    stop("Argument '", deparse(substitute(x)),
         "' must be a list.", call.=FALSE)
}
