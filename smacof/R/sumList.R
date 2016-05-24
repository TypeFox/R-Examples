#sum of distances over 

sumList <- function(x) 
{
  m <- length(x)
  z <- x[[1]]
  if (m == 1) return(z)
  for (j in 2:m) z <- z+x[[j]]
  return(z)
}
