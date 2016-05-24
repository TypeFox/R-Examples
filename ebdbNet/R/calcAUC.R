`calcAUC` <-
function(sens, cspec) 
{
    	if(is.vector(sens) != TRUE) {
        	stop("Error: ", paste(sQuote("sens"), sep = ""), " should be a vector.")
    	}
    	if(is.vector(cspec) != TRUE) {
        	stop("Error: ", paste(sQuote("cspec"), sep = ""), " should be a vector.")
    	}
	o <- order(cspec)
	data <- rbind(c(0,0),cbind(cspec[o], sens[o]), c(1,1))
	AUC <- 0
	for(i in 2:dim(data)[1]) {
		a = data[i,1] - data[(i-1),1]
		b1 = min(data[i,2], data[(i-1),2])
		b2 = max(data[i,2], data[(i-1),2]) 
		AUC = AUC + a/2*(b1+b2)
	}
	return(AUC)
}

