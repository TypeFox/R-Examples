# Description: 	A function to calculate credible intervals and make a table.  See page 130.
# Usage: 	t.ci.table(coefs,cov.mat,level=0.95,degrees=Inf,quantiles=c(0.025,0.500,0.975)) 
# Arguments: 	coefs    	vector of coefficient estimates, usually posterior means
#		cov.mat		variance-covariance matrix   
#		level		desired coverage level
#		degrees		degrees of freedom parameter for students-t distribution assumption
#		quantiles	vector of desired CDF points (quantiles) to return
# Values:	quantile.mat	matrix of quantiles


t <- function(coefs,cov.mat,level=0.95,degrees=Inf,quantiles=c(0.025,0.500,0.975))
UseMethod("t")
t.ci <- function(coefs,cov.mat,level=0.95,degrees=Inf,quantiles=c(0.025,0.500,0.975))  
{
quantile.mat <- cbind( coefs, sqrt(diag(cov.mat)),
t(qt(quantiles,degrees) %o% sqrt(diag(cov.mat)))
+ matrix(rep(coefs,length(quantiles)),
ncol=length(quantiles)) )
quantile.names <- c("Mean","Std. Error")
for (i in 1:length(quantiles))
quantile.names <- c(quantile.names,paste(quantiles[i],
"Quantile"))
dimnames(quantile.mat)[2] <- list(quantile.names)
return(title="Posterior Quantities",round(quantile.mat,4))

}


