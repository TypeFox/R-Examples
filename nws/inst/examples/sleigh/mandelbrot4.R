selectRegion <- function() {
  selected <- tryCatch(locator(2), error=function(e) NULL)
  if (is.null(selected)) return(NULL)

  x1 <- min(selected$x)
  x2 <- max(selected$x)
  y1 <- min(selected$y)
  y2 <- max(selected$y)

  if ((x2 - x1) > (y2 - y1)) {
    y2 <- y1 + (x2 - x1)
  } else {
    x2 <- x1 + (y2 - y1)
  }
  list(x1=x1, x2=x2, y1=y1, y2=y2)
}

vmandelbrot <- function(av, b0, lim) {
  mandelbrot <- function(a0) {
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
  unlist(lapply(av, mandelbrot))
}

n <- 360
m <- 100

if (dev.cur() == 1) get(getOption("device"))()

library(nws)
s <- sleigh()
accum <- function(results, indices) {
  a <- do.call(cbind, results)
  z[, indices] <<- a
  image(x, y, z, ylim=c(min(indices), max(indices)), add=TRUE, col=c(rainbow(m), '#000000'))
}
opts <- list(accumulator=accum, chunkSize=10, loadFactor=4)

x <- seq(-2.0, 0.6, length.out=n)
y <- seq(-1.3, 1.3, length.out=n)
z <- matrix(m, n, n)
image(x, y, z, main='Computing...', col='black')
eachElem(s, vmandelbrot, list(b0=y), list(av=x, lim=m), eo=opts)

xbase <- x
ybase <- y
zbase <- z

repeat {
  image(xbase, ybase, zbase, xlab='x', ylab='y',
        main='Please select a new region', col=c(rainbow(m), '#000000'))
  r <- selectRegion()
  if (is.null(r)) break
  x <- seq(r$x1, r$x2, length.out=n)
  y <- seq(r$y1, r$y2, length.out=n)
  z <- matrix(m, n, n)
  image(x, y, z, main='Computing...', col='black')
  eachElem(s, vmandelbrot, list(b0=y), list(av=x, lim=m), eo=opts)
  cat('computed:', r$x1, r$x2, ',', r$y1, r$y2, '\n')
  image(x, y, z, main='Click once to return', col=c(rainbow(m), '#000000'))
  if (is.null(tryCatch(locator(1), error=function(e) NULL))) break
}
