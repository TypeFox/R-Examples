.dist.xy <- function(x, y, npts)
{
  n <- length(x)
  
  z <- array(data = 0.0, dim = n, dimnames = NULL)

  for (i in 1:n) {
    z[i] <- sum((x <= x[i]) & (y <= y[i]))
  }
  
  z <- z / n
  
  i <- !duplicated(data.frame(x, y))

  x <- x[i] 
  y <- y[i]
  z <- z[i]  
  
  n <- length(z)
  
  if (n > npts) {
    i <- sample.int(n, npts, replace = FALSE, prob = NULL)  
  
    x <- x[i]
    y <- y[i]
    z <- z[i]
  }

  output <- list()

  output$x <- x
  output$y <- y
  output$z <- z

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .dist.xy
