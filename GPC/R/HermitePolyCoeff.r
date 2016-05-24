	HermitePolyCoeff <- function(degree) { 
		if (degree == 0){
			value <- c(1.0, 0.0, 0.0)}
		else{
			value <- c(1.0, 0.0, -(degree-1.0))}
	return(value)
	} 