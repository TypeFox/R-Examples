# a robust skewness measure described in Kim and White (2004)
rob.sk <- function(x){
  if(!is.vector(x)){
     x <- as.matrix(x)
     nobs <- dim(x)[1]; ndim <- dim(x)[2]; m <- colMeans(x)
     sk1 <- numeric(ndim)
     sk2 <- numeric(ndim)
     for(i in 1:ndim){
      # standard skewness
         x. <- x[,i]-m[i]; v <- mean(x.^2)
         std.x <- x./sqrt(v)
         sk1[i] <- mean(std.x^3)
      # robustified skewness
         quant <- quantile(x[,i], prob=seq(0.25,0.75,0.25))
         sk2[i] <- (quant[3]+quant[1]-2*quant[2])/(quant[3]-quant[1])
     }
  } else {
         x <- as.vector(x)
         nobs <- length(x); m <- mean(x); ndim <- 1
      # standard skewness
         x. <- x-m; v <- mean(x.^2)
         std.x <- x./sqrt(v)
         sk1 <- mean(std.x^3)
      # robustified skewness
         quant <- quantile(x, prob=seq(0.25,0.75,0.25))
         sk2 <- (quant[3]+quant[1]-2*quant[2])/(quant[3]-quant[1])
  }
      sk <- rbind(sk1, sk2); rownames(sk) <- c("standard", "robust"); colnames(sk) <- paste("series",1:ndim)
      sk
}
