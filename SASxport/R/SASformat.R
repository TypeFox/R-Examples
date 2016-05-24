
#' @export

SASformat <- function(x, default)
  UseMethod("SASformat")

#' @export

SASformat.default <- function(x, default=NULL)
{
  lab <- attr(x,"SASformat")
  if(is.null(lab))
    default
  else
  lab
}

#' @export

SASformat.data.frame <- function(x, default=NULL)
{
  sapply( x, SASformat)
}

#' @export

"SASformat<-" <- function(x, value)
  UseMethod("SASformat<-")

#' @export

"SASformat<-.default" <- function(x, value)
{
  attr(x,'SASformat') <- value
  x
}

#' @export

"SASformat<-.data.frame" <- function(x, value)
{
  if( ncol(x) != length(value) )
    stop("vector of formats must match number of data frame columns")
  
  for(i in 1:ncol(x))
    attr(x[[i]],'SASformat') <- value[i]
  x
}
