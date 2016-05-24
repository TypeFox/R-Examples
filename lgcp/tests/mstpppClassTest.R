library(lgcp)
library(spatstat)

xyt <- cbind(matrix(runif(200,-5,5),100,2),runif(100,0,5))
tlim <- c(0,5) 
ow <- owin(xrange=c(-5,5),yrange=c(-5,5))
data <- as.data.frame(matrix(runif(500),100,5))

PP <- ppp(xyt[,1],xyt[,2],window=ow)
ST <- stppp(list(data=xyt,tlim=tlim,window=ow))

mstppp1 <- mstppp(PP,t=xyt[,3],tlim=tlim,data=data)
mstppp2 <- mstppp(list(xyt=xyt,tlim=tlim,window=ow,data=data))
mstppp3 <- mstppp(ST,data=data)

plot(mstppp1,cols=grey(mstppp1$data[,1]/max(mstppp1$data[,1])))

pppver <- as.ppp(mstppp1)