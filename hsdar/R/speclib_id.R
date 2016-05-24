idSpeclib <- function(x)
{
if (!is.speclib(x))
  stop("Class of x must be Speclib")
ids <- x@ID

return(if (any(c(length(x@ID) != nspectra(x), anyDuplicated(ids)))) c(1:nspectra(x)) else ids)
}

"idSpeclib<-" <- function(x, value)
{
  xx <- x
  xx@ID <- value
  x <- xx
}
