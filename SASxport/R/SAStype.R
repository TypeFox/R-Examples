
#' @export

SAStype <- function(x, default)
  UseMethod("SAStype")

#' @export

SAStype.default <- function(x, default=NULL)
{
  lab <- attr(x,"SAStype")
  if(is.null(lab))
    default
  else
  lab
}

#' @export

"SAStype<-" <- function(x, value)
  UseMethod("SAStype<-")

#' @export

"SAStype<-.default" <- function(x, value)
{
  attr(x,'SAStype') <- makeSASNames(value)
  x
}
