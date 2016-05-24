head.ts <- function(x, n = 5, ...) {
  if (n > dim(as.matrix(x))[1]) {
    return(x)
  } else {
    return(window(x, start = start(x), end = start(x) + c(0, n - 1), 
      frequency=tsp(x)[3]))
  }
}

tail.ts <- function(x, n = 5, ...) {
  if (n > dim(as.matrix(x))[1]) {
    return(x)
  } else {
    return(window(x, start = end(x) - c(0, n - 1), end = end(x), 
      frequency = tsp(x)[3]))
  }
}