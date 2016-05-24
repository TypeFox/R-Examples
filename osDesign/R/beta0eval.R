beta0eval <-
function(beta0,
											betaX,
											X,
											N,
											rhoY,
											expandX)
{
	##
	if(expandX == "none") designX <- X
  ##
  if(expandX == "all")
  {
	  designX <- X[,1]
  	for(i in 2:ncol(X))
  	{
    	for(j in 1:max(X[,i]))
    	{
    		designX <- cbind(designX, as.numeric(X[,i] == j))
    	}
  	}
  }
  
  ##
  etaY  <- as.numeric(designX %*% c(beta0, betaX))
  
  ##
  value <- abs(sum(expit(etaY) * (N/sum(N))) - rhoY)
  return(value)
}
