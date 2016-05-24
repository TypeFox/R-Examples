# a robust excess kurtosis measure described in Kim and White (2004)
rob.kr <- function(x){
  if(!is.vector(x)){
     x <- as.matrix(x)
     nobs <- dim(x)[1]; ndim <- dim(x)[2]
     m <- colMeans(x)
     kr1 <- numeric(ndim)
     kr2 <- numeric(ndim)
     for(i in 1:ndim){
        x. <- x[,i]-m[i]; v <- mean(x.^2)
        std.x <- x./sqrt(v)
        kr1[i] <- mean(std.x^4)-3
        quant <- quantile(x[,i], prob=seq(0.125, 0.875, 0.125))
        kr2[i] <- (quant[7]-quant[5]+quant[3]-quant[1])/(quant[6]-quant[2])-1.23
     }
  } else {
     x <- as.vector(x)
     nobs <- length(x); m <- mean(x); ndim <- 1
     x. <- x-m; v <- mean(x.^2)
     std.x <- x./sqrt(v)
     kr1 <- mean(std.x^4)-3
     quant <- quantile(x, prob=seq(0.125, 0.875, 0.125))
     kr2 <- (quant[7]-quant[5]+quant[3]-quant[1])/(quant[6]-quant[2])-1.23
  }
  kr <- rbind(kr1, kr2); rownames(kr) <- c("standard", "robust"); colnames(kr) <- paste("series",1:ndim)
  kr
}
