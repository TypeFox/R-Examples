plot.Cuminc <- function(x, ...)
{
  if (!inherits(x, "Cuminc"))
    stop("'x' must be a 'Cuminc' object")
  
  plot(attr(x, "survfit"), ...)
}
