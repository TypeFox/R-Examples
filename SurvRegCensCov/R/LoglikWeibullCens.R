LoglikWeibullCens <- function(x, data, lowerbound, vdelta){
	if(x[1] > 0 & x[2] > 0){
		if(length(lowerbound[is.na(lowerbound)==FALSE])!=0){
			loglik <- sum(log(dweibull(data[vdelta == 1], shape=x[1], scale=x[2]))) + sum(log(pweibull(data[vdelta == 0 & is.na(lowerbound)==TRUE], shape=x[1], scale=x[2]))) + sum(log(pweibull(data[vdelta == 0 & is.na(lowerbound)==FALSE], shape=x[1], scale=x[2])-pweibull(lowerbound[vdelta == 0 & is.na(lowerbound)==FALSE], shape=x[1], scale=x[2])))
		}
		if(length(lowerbound[is.na(lowerbound)==FALSE])==0){
			loglik <- sum(log(dweibull(data[vdelta == 1], shape=x[1], scale=x[2]))) + sum(log(pweibull(data[vdelta == 0], shape=x[1], scale=x[2])))
		}
	}
     
	if(x[1] <= 0 | x[2] <= 0){ loglik <- -Inf }
     if(loglik == -Inf){loglik <- -100000}

     return(loglik)
}







#