myQQNorm <- function(x, nSim){
  n_size <- length(x)
  
  if(n_size < 5){
    stop("too small sample")
  }
  
  n_sd <- mad(x)
  n_mean <- mean(x)
  
  a <- ifelse(n_size <= 10, 3/8, 1/2)
  myseq <- seq(from = (1 - a)/(n_size + (1-a)-a), to = (n_size - a)/(n_size + (1-a)-a), length.out = min(n_size,50))
  
  n_quant <- myQQNormIntern(myseq,n_mean,n_sd,length(x),nSim)
  
  myX <- qnorm((1:n_size - a)/(n_size + (1-a)-a))
  myY <- sort(x)
  
  plot(myX, myY, main = "Normal Q-Q Plot - SIM", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles") 
  matlines(qnorm(myseq), n_quant, lty = 1,pch = 1, lwd = 3, col = "#cdd2d015")
  
  points(myX, myY)
  qqline(x)
  box()
}


setGeneric("qqnormSim", function(x, nSim = 500) 
  standardGeneric("qqnormSim"))


setMethod("qqnormSim","lm",
          function(x, nSim = 500){
            myQQNorm(resid(x), nSim)
          })

setMethod("qqnormSim","numeric",
          function(x, nSim = 500){
            myQQNorm(x, nSim)
          })
