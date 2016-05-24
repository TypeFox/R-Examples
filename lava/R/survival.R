##' @export
`survival` <- function(x,...) UseMethod("survival")
##' @export
"survival<-" <- function(x,...,value) UseMethod("survival<-")

##' @export
"survival<-.lvm" <- function(x,...,value) {
  survival(x, value, ...)
}

##' @export
`survival.lvm` <-
function(x,var=NULL, ...) {
  if (is.null(var)) {
    survidx <- unlist(x$attributes$survival)
    if (length(survidx)>0)
      return(names(survidx)[survidx])
    else
      NULL
  }
  x$attributes$survival[var] <- TRUE
  x$attributes$normal[var] <- TRUE
  return(x)
}
