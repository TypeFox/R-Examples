`print.dlmap` <- 
function(x, ...)
{
  cat(" This is an object of class \"dlmap\".\n")
  cat(" It is too complex to print, so we provide a summary of its input data and final results.\n")
  summary(x)
}
