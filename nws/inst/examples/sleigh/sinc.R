x <- seq(-10, 10, length=200)
f <- function(x, y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }
library(nws)
s <- sleigh()
z <- unlist(eachElem(s, f, expand.grid(x=x, y=x), eo=list(chunkSize=200)))
dim(z) <- c(length(x), length(x))
persp(x, x, z, ylab='y', theta=30, phi=30, expand=0.5, col="lightblue")
