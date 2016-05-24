greplace <-
function(data,pattern,replacement,...){
	if(length(replacement)==1){
		n <- length(pattern)
		for(i in 1:n)data[data==pattern[i]] <- replacement
	}
	else{
		if(length(replacement)!= length(pattern))stop("The mumber of pattern is not equal to the number of replacement!")
		else {
			inter <- intersect(pattern,replacement)
			if(length(inter)==0){
				n <- length(pattern)
				for(i in 1:n)data[data==pattern[i]] <- replacement[i]
			}
			if(length(inter)>0){
				dat <- data
				n <- length(pattern)
				for(i in 1:n)dat[data==pattern[i]] <- replacement[i]
				data <- dat
			}
		}
	}
	return(data)	
  }
