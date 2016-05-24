jOptWidthKernel <-
function(s){
# s is data 
 n <- length(s)
c1 <- 1
c2 <- 1/2/sqrt(pi)
#tinhs c3:
c30 <- 1
 h <- c1^(-2/5)*c2^(1/5)*c30^(-1/5)/n^(1/5)
dc <- 1
numloop <- 0
while (dc > 0.01 & numloop < 100) {
c3 <- jSimsonKernel(1000,-10,10,s,h,jGradKernelFitPdf2)
h <- c1^(-2/5)*c2^(1/5)*c3^(-1/5)/n^(1/5)
dc <- c3 - c30
c30 <- c3
numloop <- numloop + 1
}
return(h)
}
