InfoMatrix.beta <- function(x,weight){

#create vector of ones
   if(is.matrix(x)==TRUE || is.data.frame(x) == TRUE) {onevec <- rep(1,length(x[,1]))} else onevec <- rep(1,length(x)) 

#create real X matrix, namely where 1st col is vector of 1s to incorporate beta0 effect
   newx <- cbind(onevec,x) 
   xmat <- as.matrix(newx)

# Create the I(beta) block to the information matrix
   W <- diag(weight)
   I <- t(xmat) %*% W %*% xmat

return(I)
}

