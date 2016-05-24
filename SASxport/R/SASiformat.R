
#' @export

SASiformat <- function(x, default)
  UseMethod("SASiformat")

#' @export

SASiformat.default <- function(x, default=NULL)
{
  lab <- attr(x,"SASiformat")
  if(is.null(lab))
    default
  else
  lab
}

#' @export

"SASiformat<-" <- function(x, value)
  UseMethod("SASiformat<-")

#' @export

"SASiformat<-.default" <- function(x, value)
{
  attr(x,'SASiformat') <- value
  x
}
