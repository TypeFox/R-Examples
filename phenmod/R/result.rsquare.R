result.rsquare <- function(values, type="cod"){
	values.observed <- values$doy.observed
	values.calculated <- values$doy.model
	check <- which((values.observed > -9999) & (!is.na(values.observed)) & 
				(values.calculated > -9999) & (!is.na(values.calculated)))

	values.observed <- values.observed[check]
	values.calculated <- values.calculated[check]

	if (type=="cod"){
		total.sum <- sum( (values.observed - mean(values.observed))^2 )
		residual.sum <- sum( (values.observed - values.calculated)^2 )

		rsquare <- 1 - (residual.sum / total.sum )
	} else {
		if (type=="pearson"){
			rsquare <- (sum( (values.observed - mean(values.observed))*
					(values.calculated - mean(values.calculated)) ) /
				(sqrt(sum( (values.observed - mean(values.observed))^2 ))*
					sqrt(sum( (values.calculated - mean(values.calculated))^2 ))))^2
			
		} else {
			rsquare <- NA
			cat("Wrong type for rsquare-calculation!\n",sep="")
		}
	}

	return(rsquare)
}
