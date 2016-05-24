library(minerva)
x <- 0:200 / 200
y <- sin(10 * pi * x) + x
mine(x,y, n.cores=1)


