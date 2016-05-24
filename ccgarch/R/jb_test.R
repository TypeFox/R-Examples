# Lomnicki(1961)-Jarque-Bera (1987) test of normality
jb.test <- function(x){
  if(!is.vector(x)){
      nobs <- dim(x)[1]; ndim <- dim(x)[2]; m <- colMeans(x)
      jb <- numeric(ndim)  
         for(i in 1:ndim){
            x. <- x[,i]-m[i]; v <- mean(x.^2)
            std.x <- x./sqrt(v)
            sk <- mean(std.x^3); kr <- mean(std.x^4)-3
            jb[i] <- nobs/6*(sk^2+0.25*kr^2)
         }
  } else {
      nobs <- length(x); ndim <- 1
      x. <- x-mean(x); v <- mean(x.^2)
      std.x <- x./sqrt(v)
      sk <- mean(std.x^3); kr <- mean(std.x^4)-3
      jb <- nobs/6*(sk^2+0.25*kr^2)
  }
      p.val <- pchisq(jb,2,lower.tail=FALSE)
      out <- rbind(jb, p.val)
      rownames(out) <- c("test stat", "p-value")
      colnames(out) <- paste("series",1:ndim)
      out
}
