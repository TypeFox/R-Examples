informationMatrixARTFIMA <- function(d, lambda){
  v11 <- 4*pi*exp(-2*lambda)*(1+ 0.5*exp(-2*lambda)) #eqn (43)
  v22 <- (d^2)*exp(-2*lambda)/(1-exp(-2*lambda)) #eqn (48)
  v12 <- d*log(1-exp(-2*lambda)) #eqn (54)
  matrix(c(v11,v12,v12,v22), ncol=2, 
         dimnames=list(c("lambda","d"),c("lambda","d")))
}

seARTFIMA <- function(d, lambda, n){
  v <- solve(informationMatrixARTFIMA(d, lambda)*n)
  sqrt(diag(v))
}