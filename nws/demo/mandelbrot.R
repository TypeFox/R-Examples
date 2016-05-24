vmandelbrot <- function(av, b0, lim) {
  mandelbrot <- function(a0) {
    a <- a0; b <- b0; a2 <- a * a; b2 <- b * b
    n <- 0
    while (a2 + b2 < 4 && n < lim) {
      b <- 2 * a * b + b0; a <- a2 - b2 + a0
      a2 <- a * a; b2 <- b * b
      n <- n + 1
    }
    n
  }
  unlist(lapply(av, mandelbrot))
}

x <- seq(0.33, 0.36, length.out=240)
y <- seq(-0.085, -0.055, length.out=240)
m <- 100

library(nws); s <- sleigh()
z <- matrix(-1, length(x), length(y))
breaks <- c(-1:62 + 0.5, m - 0.5, 1000000)
col <- c(rainbow(64), '#000000')
image(x, y, z, main='Processing with 3 local workers...', breaks=breaks, col=col)

accum <- function(results, indices) {
  a <- do.call(cbind, results)
  z[, indices] <<- a
  image(x, y, z, ylim=c(min(indices), max(indices)), add=TRUE, breaks=breaks, col=col)
}
opts <- list(accumulator=accum, chunkSize=10, loadFactor=4)
eachElem(s, vmandelbrot, list(b0=y), list(av=x, lim=m), eo=opts)
image(x, y, z, main='Mandelbrot Demo Complete', breaks=breaks, col=col)
stopSleigh(s)
