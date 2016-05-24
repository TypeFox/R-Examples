Reciprocal <- function(x){
	
	n = length(x)
	temp = c()
	for(i in 1: n){
		if(x[i] == 0) temp[i] = 0
		else temp[i] = 1/x[i]
	}
	return(temp)
}