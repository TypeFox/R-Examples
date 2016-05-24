
predict_update_km_parallel <- function(newXmean,newXvar,newXvalue, Sigma.r,
		newdata.oldmean,newdata.oldsd,kn){
			
	chol.Sigma.r <- NULL
	chol.Sigma.r <- try(chol(Sigma.r),TRUE)
	if(!is.numeric(chol.Sigma.r)) return(list(error=TRUE))
	
	lambda_nplus.r <- kn %*% chol2inv(chol.Sigma.r)
	newXdiff <- matrix(newXvalue - newXmean,ncol=1)
	
	predict_mean <- newdata.oldmean + lambda_nplus.r %*% newXdiff 
	predict_var <- pmax(0,newdata.oldsd*newdata.oldsd - rowSums(lambda_nplus.r * kn))
		
	predict_sd <- sqrt(predict_var)
	
	output.list <- list(mean = predict_mean,sd=predict_sd, lambda=lambda_nplus.r)
	return(output.list)
}

