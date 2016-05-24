predict_update_km <- function(newXmean,newXvar,newXvalue, 
								newdata.oldmean,newdata.oldsd,kn){
	
	lambda_nplus1 <- kn/newXvar
	
	predict_mean <- newdata.oldmean + lambda_nplus1 * ( newXvalue - newXmean  )
	predict_var <- pmax(0,newdata.oldsd*newdata.oldsd - lambda_nplus1*lambda_nplus1 * newXvar)
	predict_sd <- sqrt(predict_var)
	
	output.list <- list(mean = predict_mean,sd=predict_sd, lambda=lambda_nplus1)
	return(output.list)
}

