getChunker <- function(numTasks, numChunks) {
  chunkSize <- numTasks %/% numChunks
  rem <- numTasks - (chunkSize * numChunks)
  f <- if (rem > 0)
    factor(c(rep(1:rem, each=chunkSize+1),
             rep((rem+1):numChunks, each=chunkSize)))
  else
    factor(rep(1:numChunks, each=chunkSize))

  function(tasks) {
    if (length(tasks) != numTasks) stop('incorrect number of tasks')
    split(tasks, f)
  }
}

x <- seq(-10, 10, length.out=200)
d <- expand.grid(x=x, y=x)
chunker <- getChunker(40000, 200)
xlist <- chunker(d$x)
ylist <- chunker(d$y)

f <- function(x, y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }
library(nws)
s <- sleigh()
z <- unlist(eachElem(s, f, elementArgs=list(x=xlist, y=ylist)))
dim(z) <- c(length(x), length(x))
persp(x, x, z, ylab='y', theta=30, phi=30, expand=0.5, col="lightblue")
