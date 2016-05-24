library("pdc")
 
fuzz <- 1e-6
if (abs(sum(pdclust(matrix(sin(1:1000)^2,ncol=4))$D)-0.150942) > fuzz) stop("test failed!") 
