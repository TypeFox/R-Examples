print.Cuminc <- function(x, ...)
{
  if (!inherits(x, "Cuminc"))
    stop("'x' must be a 'Cuminc' object")
  
  print(as.data.frame(x))
}
