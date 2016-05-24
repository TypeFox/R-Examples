simEcm <-
function(n, s){
	
	dd <- sample(n, s, replace=TRUE)
	
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
