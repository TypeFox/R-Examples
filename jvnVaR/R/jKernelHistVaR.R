jKernelHistVaR <-
function(s, alpha){
hh <- jOptWidthKernel(s)
h <- 0.01
x <- min(s)
cdf <- jSimsonKernel(1000,-10,x,s,hh,jKernelFitPdf)
while (cdf < alpha){
x <- x+h
cdf <- jSimsonKernel(1000,-10,x,s,hh,jKernelFitPdf)
}

return(-x)
}
