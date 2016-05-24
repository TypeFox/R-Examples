

finiterange <- function (x) ## NOEXPORT
{
  x <- x[! is.na(x)]
  x <- x[is.finite(x)]
  if (length(x)==0) {
    return(c(NA, NA))
  }
  range(x)  
}

