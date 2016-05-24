f1 <- function(x, y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }
f2 <- function(x, y) { 1 - exp(-1/(x^2+y^2)) }
f3 <- function(x, y) { sin(x) * cos(y) }
f4 <- function(x, y) { 0.1 * x * sin(2*y) }

library(nws)
s <- sleigh()

if (dev.cur() == 1) get(getOption("device"))(height=10, width=10)
par(mfcol=c(2,2))
colors <- rainbow(4)

# element by element method
x <- seq(-10, 10, length=60)
t <- system.time(z <- eachElem(s, f1, expand.grid(x=x, y=x)))
z <- unlist(z)
dim(z) <- c(length(x), length(x))
persp(x, x, z, ylab='y', theta=30, phi=30, expand=0.5, col=colors[1])
title(main='Element by element', sub=paste(t[3], 'sec'))

# element by row method
x <- seq(-5, 5, length=60)
t <- system.time(z <- eachElem(s, f2, x, x))
z <- do.call(rbind, z)
persp(x, x, z, ylab='y', theta=30, phi=30, expand=0.5, col=colors[2])
title(main='Element by row', sub=paste(t[3], 'sec'))

# element by element with task chunking
x <- seq(-pi/2, pi/2, length=60)
chunkSize <- 60
opts <- list(chunkSize=chunkSize)
t <- system.time(z <- eachElem(s, f3, expand.grid(x=x, y=x), eo=opts))
z <- unlist(z)
dim(z) <- c(length(x), length(x))
persp(x, x, z, ylab='y', theta=30, phi=30, expand=0.5, col=colors[3])
title(main=paste('Element by element with chunkSize:', chunkSize), sub=paste(t[3], 'sec'))

# element by row method with task chunking
x <- seq(0, 16, length=60)
chunkSize <- 5
opts <- list(chunkSize=chunkSize)
t <- system.time(z <- eachElem(s, f4, x, x, eo=opts))
z <- do.call(rbind, z)
persp(x, x, z, ylab='y', theta=30, phi=30, expand=0.5, col=colors[4])
title(main=paste('Element by row with chunkSize:', chunkSize), sub=paste(t[3], 'sec'))
