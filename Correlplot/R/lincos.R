lincos <- function(x) {
   if(is.vector(x)) {
      y <- rep(NULL,length(x))
      x <- abs(x)
      y[x < pi]  <- -2*x[x < pi]/pi +1
      y[x >= pi] <- 2*x[x >= pi]/pi - 3
   }
   if(is.matrix(x)) {
      y <- matrix(rep(0,ncol(x)*nrow(x)),ncol=ncol(x))
      x <- abs(x)      
      y[x < pi]  <- -2*x[x < pi]/pi +1
      y[x >= pi] <- 2*x[x >= pi]/pi - 3
   }
   return(y)
}
