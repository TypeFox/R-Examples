.densParzenWindow.xy <- function(x, y, hx, hy, npts)
{
  output <- .C("RdensParzenWindowXY",
    n = as.integer(length(x)),
    x = as.double(x),
    y = as.double(y),
    z = double(length(x)),
    hx = as.double(hx),
    hy = as.double(hy),
    error = integer(1),
    PACKAGE = "rebmix")

  if (output$error == 1) {
    stop("in densParzenWindow.xy!", call. = FALSE); return(NA)
  }

  i <- !duplicated(data.frame(output$x, output$y))

  output$x <- output$x[i] 
  output$y <- output$y[i]
  output$z <- output$z[i]  
  
  n <- length(output$z)
  
  if (n > npts) {
    i <- sample.int(n, npts, replace = FALSE, prob = NULL)  
  
    output$x <- output$x[i]
    output$y <- output$y[i]
    output$z <- output$z[i]
  }

  rm(list = ls()[!(ls() %in% c("output"))])

  return(output)
} ## .densParzenWindow.xy
