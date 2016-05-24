LoglikGammaCens <- function(x, data, lowerbound, vdelta){
	if(x[1] > 0 & x[2] > 0){
		if(length(lowerbound[is.na(lowerbound) == FALSE]) != 0){
			loglik <- sum(log(dgamma(data[vdelta == 1], shape = x[1], rate = x[2]))) + sum(log(pgamma(data[vdelta == 0 & is.na(lowerbound)==TRUE], shape=x[1], rate=x[2]))) + sum(log(pgamma(data[vdelta == 0 & is.na(lowerbound)==FALSE], shape = x[1], rate = x[2]) - pgamma(lowerbound[vdelta == 0 & is.na(lowerbound) == FALSE], shape = x[1], rate = x[2])))
		}
		if(length(lowerbound[is.na(lowerbound) == FALSE]) == 0){
			loglik <- sum(log(dgamma(data[vdelta == 1], shape = x[1], rate = x[2]))) + sum(log(pgamma(data[vdelta == 0], shape = x[1], rate = x[2])))
		}
	}
	if(x[1] <= 0 | x[2] <= 0){loglik <- -Inf}
     if(loglik == -Inf){loglik <- -100000}
     
return(loglik)
}







#