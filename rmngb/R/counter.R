makeCounter <- function(start = 1) {
  i <- start
  function(f = function(x) x, ...)
    (i <<- f(i, ...))
}

increment <- function(x, step = 1)
  x + step

decrement <- function(x, step = 1)
  x - step
