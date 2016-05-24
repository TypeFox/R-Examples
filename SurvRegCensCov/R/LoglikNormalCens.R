LoglikNormalCens <- function(x, data, lowerbound, vdelta){
	if(x[2] > 0){
		if(length(lowerbound[is.na(lowerbound)==FALSE])!=0){
			loglik <- sum(log(dnorm(data[vdelta == 1], mean=x[1], sd=x[2]))) + sum(log(pnorm(data[vdelta == 0 & is.na(lowerbound)==TRUE], mean=x[1], sd=x[2]))) + sum(log(pnorm(data[vdelta == 0 & is.na(lowerbound)==FALSE], mean=x[1], sd=x[2])-pnorm(lowerbound[vdelta == 0 & is.na(lowerbound)==FALSE], mean=x[1], sd=x[2])))
		}
		if(length(lowerbound[is.na(lowerbound)==FALSE])==0){
			loglik <- sum(log(dnorm(data[vdelta == 1], mean=x[1], sd=x[2]))) + sum(log(pnorm(data[vdelta == 0], mean=x[1], sd=x[2])))
		}
	}
	
     if(x[2] <= 0){ loglik <- -Inf }
     if(loglik == -Inf){loglik <- -100000}

     return(loglik)
}







#