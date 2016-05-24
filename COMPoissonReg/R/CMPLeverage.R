CMPLeverage <- function(x,y,betahat,nuhat,max){

#  create vector of ones
   if(is.matrix(x)==TRUE || is.data.frame(x) == TRUE) {onevec <- rep(1,length(x[,1]))} else onevec <- rep(1,length(x)) 

#  create real X matrix, namely where 1st col is vector of 1s to incorporate beta0 effect
   newx <- cbind(onevec,x) 
   xmat <- as.matrix(newx)

# 1) to code the W matrix  (diagonal matrix with Var(Y_i) )

   W <- diag(weights(x,betahat,nuhat,max))

#    and X matrix (in Appendix)

   E.y <- computez.prodj(xmat,betahat,nuhat,max)/computez(xmat,betahat,nuhat,max)
   E.logfacty <- computez.prodlogj(xmat,betahat,nuhat,max)/computez(xmat,betahat,nuhat,max)
   extravec <- (-log(factorial(y)) + E.logfacty)/(y - E.y)
   curlyX.mat <- cbind(xmat,extravec)

# 2) to compute H using eq (12)  on p. 11

   H1 <- t(curlyX.mat) %*% sqrt(W)
   H2 <- solve(t(curlyX.mat) %*% W %*% curlyX.mat)
   H <- t(H1) %*% H2 %*% H1
   diagH <- diag(H)
return(diagH)
}

