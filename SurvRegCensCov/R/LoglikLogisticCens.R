LoglikLogisticCens <- function(x, data, lowerbound, vdelta){
	if(x[2] > 0){
		if(length(lowerbound[is.na(lowerbound)==FALSE])!=0){
			loglik <- sum(log(dlogis(data[vdelta == 1], location=x[1], scale=x[2]))) + sum(log(plogis(data[vdelta == 0 & is.na(lowerbound)==TRUE], location=x[1], scale=x[2]))) + sum(log(plogis(data[vdelta == 0 & is.na(lowerbound)==FALSE], location=x[1], scale=x[2])-plogis(lowerbound[vdelta == 0 & is.na(lowerbound)==FALSE], location=x[1], scale=x[2])))
		}
		if(length(lowerbound[is.na(lowerbound)==FALSE])==0){
			loglik <- sum(log(dlogis(data[vdelta == 1], location=x[1], scale=x[2]))) + sum(log(plogis(data[vdelta == 0], location=x[1], scale=x[2])))
		}
	}
	
     if(x[2] <= 0){ loglik <- -Inf }
     if(loglik == -Inf){loglik <- -100000}
     
     return(loglik)
}







#