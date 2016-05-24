constantCMPfitsandresids <- function(betahat,nuhat,x,y=0){

#create vector of ones
   if(is.matrix(x)==TRUE || is.data.frame(x) == TRUE) {onevec <- rep(1,length(x[,1]))} else onevec <- rep(1,length(x)) 

#create real X matrix, namely where 1st col is vector of 1s to incorporate beta0 effect
   newx <- cbind(onevec,x) 
   xmat <- as.matrix(newx)

#Determine estimated lambdahat, fit, and residuals
   lambdahat <- exp(xmat %*% betahat)
   fit <- lambdahat^(1/nuhat) - ((nuhat -1)/(2*nuhat))

   resid <- y - fit

return(list(fit=fit,resid=resid))
}

