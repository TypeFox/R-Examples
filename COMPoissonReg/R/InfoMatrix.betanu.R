InfoMatrix.betanu <- function(x,beta,nu,max){

# create vector of ones
   if(is.matrix(x)==TRUE || is.data.frame(x) == TRUE) {onevec <- rep(1,length(x[,1]))} else onevec <- rep(1,length(x)) 

# create real X matrix, namely where 1st col is vector of 1s to incorporate beta0 effect
   newx <- cbind(onevec,x) 
   xmat <- as.matrix(newx)

# Compute the components needed to determine the beta,nu block of the information matrix
  ans1 <- computez.prodjlogj(xmat,beta,nu,max)
  ans2 <- computez.prodlogj(xmat,beta,nu,max)
  ans3 <- computez.prodj(xmat,beta,nu,max)
  ans4 <- computez(xmat,beta,nu,max)

  ans <- (ans1/ans4) - ((ans2/ans4) * (ans3/ans4))

# Compute the beta,nu block of the information matrix
  I <- t(xmat) %*% ans

return(I)
}

