.dist.x <- function(x, npts)
{
  n <- length(x)
  
  y <- array(data = 0.0, dim = n, dimnames = NULL)

  for (i in 1:n) {
    y[i] <- sum(x <= x[i])
  }
  
  y <- y / n
  
  i <- !duplicated(x) 

  x <- x[i] 
  y <- y[i]
  
  n <- length(y)
  
  if (n > npts) {
    i <- sample.int(n, npts, replace = FALSE, prob = NULL)  
  
    x <- x[i]
    y <- y[i]
  }

  output <- list()

  output$x <- x
  output$y <- y

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .dist.x
