.densKNearestNeighbour.x <- function(x, k, hx, npts)
{
  output <- .C("RdensKNearestNeighbourX",
    n = as.integer(length(x)),
    x = as.double(x),
    y = double(length(x)),
    k = as.integer(k),
    hx = as.double(hx),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in densKNearestNeighbour.x!", call. = FALSE); return(NA)
  }

  i <- !duplicated(output$x) 

  output$x <- output$x[i] 
  output$y <- output$y[i]
  
  n <- length(output$y)
  
  if (n > npts) {
    i <- sample.int(n, npts, replace = FALSE, prob = NULL)  
  
    output$x <- output$x[i]
    output$y <- output$y[i]
  }
  
  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densKNearestNeighbour.x
