print.mpqtl <-
function(x, ...)
{
  cat(" This is an object of class \"mpqtl\".\n")
  cat(" It is too complex to print, so we provide the following summary.\n")
  print(summary(x))
}

