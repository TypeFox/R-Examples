normalPdf <-
function(x, mean, sd){
	a = 1/(sd*sqrt(2*pi))
	b = exp(-((x-mean)^2)/(2*sd^2))
	return(a*b)
	}
