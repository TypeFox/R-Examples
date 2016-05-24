LRT <- function(x,y,betahat0,betahat,nuhat,max){
# This function computes the -2logLRT value and associated p-value for significance

#create vector of ones
  if(is.matrix(x)==TRUE || is.data.frame(x) == TRUE) {onevec <- rep(1,length(x[,1]))} else onevec <- rep(1,length(x)) 

#create real X matrix, namely where 1st col is vector of 1s to incorporate beta0 effect
  newx <- cbind(onevec,x) 
  xmat <- as.matrix(newx)

# Compute the test statistic
  teststat <- -2*((t(y)%*%(xmat %*% betahat0)) - sum(log(factorial(y))) - sum(exp(xmat %*% betahat0)) - 
              ((t(y)%*%(xmat %*% betahat)) - (nuhat*sum(log(factorial(y)))) - 
              sum(log(computez(xmat,betahat,nuhat,max)))))

# Determine the associated p-value
  pvalue <- pchisq(teststat,df=1,lower.tail=FALSE)

return(list(teststat=teststat,pvalue=pvalue))
}

