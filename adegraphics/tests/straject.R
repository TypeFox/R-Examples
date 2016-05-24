library(adegraphics)
pdf("straject.pdf")

rw <- function(a){
  x <- 0
  for(i in 1:49)
    x <- c(x, x[length(x)] + runif(1, -1, 1))
  x
}

x <- unlist(lapply(1:5, rw))
y <- unlist(lapply(1:5, rw))
z <- gl(5, 50)
g1 <- s.traject(data.frame(x, y), z, ppoints.pch = 19:23, plines.col = rainbow(5))

x <- unlist(lapply(1:2, rw))
y <- unlist(lapply(1:2, rw))
z <- gl(2, 50)
g2 <- s.traject(data.frame(x, y), z, ppoints.pch = 21:20, plines.col = 1:2)
