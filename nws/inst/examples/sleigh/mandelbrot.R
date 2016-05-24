mandelbrot <- function(a0, b0, lim) {
  a <- a0; b <- b0
  a2 <- a * a
  b2 <- b * b
  n <- 0
  while (a2 + b2 < 4 && n < lim) {
    b <- 2 * a * b + b0
    a <- a2 - b2 + a0
    a2 <- a * a     
    b2 <- b * b
    n <- n + 1
  }
  n
}

ea <- function() {
  ix <- iy <- 1
  nx <- length(x)
  ny <- length(y)

  function() {
    if (iy > ny) {
      stop()
    } else {
      r <- list(x[ix], y[iy], m)
      ix <<- ix + 1
      if (ix > nx) {
        ix <<- 1
        iy <<- iy + 1
      }
      r
    }
  }
}

x <- seq(-2.0, 0.6, length.out=240)
y <- seq(-1.3, 1.3, length.out=240)
m <- 100

library(nws)
s <- sleigh()
opts <- list(elementFunc=ea(), chunkSize=length(x)*10, loadFactor=4)
z <- unlist(eachElem(s, mandelbrot, eo=opts))
dim(z) <- c(length(x), length(y))

if (dev.cur() == 1) get(getOption("device"))()
image(x, y, z, col=c(rainbow(m), '#000000'))
