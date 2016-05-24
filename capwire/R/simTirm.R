simTirm <-
function(na, nb, alpha, s){
	
	v <- c(1:(na+nb))
	dd <- sample(v, s, prob=c(rep(1,na), rep(alpha, nb)), replace=TRUE)
	
	data <- data.frame()
	ind <- unique(dd)
	for (i in 1:length(ind)){
		x <- length(dd[dd == ind[i]])
		y <- c(ind[i], x)
		data <- rbind(data, y)
	}

	data <- indToClass(data)
	
	return(data)	
}
