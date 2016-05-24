#relative influence function, cf: Friedman (2001, Li and Luan (2005)
# obtained from help(predict.smooth.spline)

     ## "Proof" that the derivatives are okay, by comparing with approximation
     diff.quot <- function(x,y) {
       ## Difference quotient (central differences where available)
       n <- length(x); i1 <- 1:2; i2 <- (n-1):n
       c(diff(y[i1]) / diff(x[i1]), (y[-i1] - y[-i2]) / (x[-i1] - x[-i2]),
         diff(y[i2]) / diff(x[i2]))
     }
#x: input design vector, %row is observation and col is covariates
#y: predicted value
rif <- function(x,y){
# res <- rep(NA,ncol(x))
# for (i in 1:ncol(x))
#  res[i] <- mean((diff.quot(x[,i], y)^2))*var(x[,i])
  res <- mean((diff.quot(x, y)^2))*var(x)
 res <- sqrt(res)
 return(res)
}   

