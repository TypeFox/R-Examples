bandnames <- function(x)
{
  if (is.null(attr(x, "bandnames")))
  {
    return(paste("B", wavelength(x), sep = "_"))
  } else {
    return(attr(x, "bandnames"))
  }
}


"bandnames<-" <- function(x, value)
{
  xx <- x
  attr(xx, "bandnames") <- value
  x <- xx
}