
## Compute AUC

AUC.calc <- function(at, values){
	AUC <- 0
	unique <- unique(at)
	values2 <- c()
	for(i in 1:length(unique)){
		values2[i] <- median(values[which(at==unique[i])])
	}
	if(unique[1]> unique[length(unique)]){
		unique <- unique[length(unique):1]
		values2<- values2[length(values2):1]
	}
	x <- c(0,unique,1)	
	y <- c(0,values2,1)
	for(i in 1:(length(y)-1)){
		AUC <- (x[i+1]-x[i])*((1/2)*(y[i+1]+y[i])) + AUC
	}
	return (AUC)
}

