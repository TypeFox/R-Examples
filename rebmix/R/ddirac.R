ddirac <- function(x, location)
{
  f <- array(data = 0.0, dim = length(x), dimnames = NULL)

  for (i in 1:length(x)) {
    if (x[i] == location) {
      f[i] = 1.0
    }
    else {
      f[i] = 0.0
    }
  }

  rm(list = ls()[!(ls() %in% c("f"))])
 
  return(f)
} ## ddirac
