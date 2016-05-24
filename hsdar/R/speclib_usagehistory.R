usagehistory <- function (x)
  return(x@usagehistory)

"usagehistory<-" <- function(x, value)
{
  xx <- x
  if (is.null(value))
  {
    xx@usagehistory <- character()
  } else {
    xx@usagehistory <- if (is.null(xx@usagehistory)) value else c(xx@usagehistory, value)
  }
  x <- xx
}