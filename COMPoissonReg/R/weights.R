weights <- function(x,beta,nu,max){

#create vector of ones
  if(is.matrix(x)==TRUE || is.data.frame(x) == TRUE) {onevec <- rep(1,length(x[,1]))} else onevec <- rep(1,length(x)) 

#create real X matrix, namely where 1st col is vector of 1s to incorporate beta0 effect
  newx <- cbind(onevec,x) 
  xmat <- as.matrix(newx)

# Compute the parts that comprise the weight functions
  w1 <- computez.prodj2(xmat,beta,nu,max)
  w2 <- computez.prodj(xmat,beta,nu,max)
  w3 <- computez(xmat,beta,nu,max)

  Ey2 <- w1/w3
  E2y <- (w2/w3)^2

# Determine the weights
  weight <- Ey2 - E2y

return(weight)
}

