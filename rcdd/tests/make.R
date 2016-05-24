
library(rcdd)

x <- seq(10, 90, 10)
x <- x[x != 50]
y <- as.numeric(x > 50)
m <- cbind(1, x)

tanv <- m
tanv[y == 1, ] <- (-tanv[y == 1, ])
tanv

makeV(points = tanv)
makeV(rays = tanv)
makeV(lines = tanv)

