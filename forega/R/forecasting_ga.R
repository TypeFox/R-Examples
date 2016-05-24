require("Rcpp")
require("robfilter")
require("forecast")
#sourceCpp("forecasting_ga.cpp")

ForecastArima = function(genes, Anchestors, MinimumForecastLength = 5, forecastprob, robust.filter.width=5){
    operator.applied <- FALSE
    n <- dim(Anchestors)[1]
    p <- dim(Anchestors)[2]
    result <- genes
    if(n < MinimumForecastLength){
        return(result)
    }
    for (i in 1:p){
        if(runif(1) < forecastprob){
			A <- Anchestors[,i]
			A <- A[!is.nan(A)]
			A <- A[!is.na(A)]
			if(length(A)<=MinimumForecastLength){
				return(result)
			}
            tseries1 <- unlist(rm.filter(A, robust.filter.width)$level)
			tseries1 <- tseries1[!is.nan(tseries1)]
			tseries1 <- tseries1[!is.na(tseries1)]
			#print(tseries1)
		
			result = tryCatch({
	        	model1 <- auto.arima(tseries1)
            	f1 <- forecast.Arima(model1,1)$mean
           		if(!is.na(f1) || !is.nan(f1)){
					result[i] <- f1
				}else{
					result[i] <- genes[i]
				}
			}, warning = function(w) {
			    #cat("W")
			}, error = function(e) {
    			result[i] <- genes[i]
			}, finally = {
			    #cat("E")
			})
			
        #cat(".")
        }
    }
	if(any(!is.na(result)) || any(!is.nan(result))){
		result <- genes
	}
    return(result);
}

