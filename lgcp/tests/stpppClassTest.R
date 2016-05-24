library(lgcp)
library(spatstat)

xyt <- cbind(matrix(runif(200,-5,5),100,2),runif(100,0,5))
tlim <- c(0,5) 
ow <- owin(xrange=c(-5,5),yrange=c(-5,5))
win <- t(sapply((0:4)*(2*pi/5),function(theta){return(5*c(cos(theta),sin(theta)))}))

stppp1 <- stppp(list(data=xyt,tlim=tlim,window=ow))

stppp2 <- stppp(list(data=xyt,tlim=tlim,window=owin(poly=win)))
