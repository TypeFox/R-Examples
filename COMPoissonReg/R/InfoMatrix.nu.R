InfoMatrix.nu <- function(x,beta,nu,max){

#create vector of ones
   if(is.matrix(x)==TRUE || is.data.frame(x) == TRUE) {onevec <- rep(1,length(x[,1]))} else onevec <- rep(1,length(x)) 

#create real X matrix, namely where 1st col is vector of 1s to incorporate beta0 effect
   newx <- cbind(onevec,x) 
   xmat <- as.matrix(newx)

# Compute the components needed to determine the information component due to nu
  ans1 <- computez.prodlogj2(xmat,beta,nu,max)
  ans2 <- computez.prodlogj(xmat,beta,nu,max)
  ans3 <- computez(xmat,beta,nu,max)

# Compute the information associated with nu
  ans <- ans1/ans3 - ((ans2/ans3)^2)

return(sum(ans))
}

