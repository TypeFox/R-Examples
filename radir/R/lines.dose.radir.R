lines.dose.radir <-
function(x, ...)
{
  if (class(x)!="dose.radir") stop("Wrong object")
  
  lines(x[[2]], x[[1]], ...)
}
