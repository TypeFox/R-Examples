
.DF.N.get <- function(data, tree){   
	temp <- sapply(tree,length)

	DF <-sum(temp)-length(temp)
	N <- rep(NA,length(temp))
	   
	for (i in 1:length(temp)){
		if (i==1) N[1] <- sum(data[1:temp[1]])
		else N[i] <- sum(data[(sum(temp[1:(i-1)])+1):sum(temp[1:i])])
	}
	return(list(DF,N))
}
