pdirac <- function(q, location, lower.tail = TRUE, log.p = FALSE)
{
  f <- array(data = 0.0, dim = length(q), dimnames = NULL)

  if (lower.tail == TRUE) {
    for (i in 1:length(q)) {
      if (q[i] >= location) {
        f[i] = 1.0
      }
      else {
        f[i] = 0.0
      }
    }
  }
  else {
    for (i in 1:length(q)) {
      if (q[i] > location) {
        f[i] = 1.0
      }
      else {
        f[i] = 0.0
      }
    }  
  }
  
  if (log.p == TRUE) {
    f <- log(f)
  }

  rm(list = ls()[!(ls() %in% c("f"))])
 
  return(f)
} ## pdirac
