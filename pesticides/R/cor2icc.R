cor2icc <-
function(x, n, type=c('within', 'not')){
	type <- tolower(type)
	if(class(x) %in% c('matrix', 'data.frame')){
		if(dim(x)[2] == 2){
			n <- dim(x)[1]
			if(n < 2){
				stop("Too few rows in x")
			}
			r <- cor(x[,1],x[,2])
		} else {
			stop("Improper number of columns in x")
		}
	} else if(class(x) != 'numeric'){
		stop("x must be a correlation or a data matrix")
	} else {
		r <- x
	}
  #	Safely assume x is a correlation
	if(type[1] %in% c('within', 'a', 'apple', 'pepper')){
		return((n*r^2-1)/(n-1))
	} else if(type[1] %in% c('not', 'peach', 'pear')){
		temp <- (n-1)*r*r + sqrt(r^4*(n-1)^2+4^r*r*n)
		return(temp/2/n)
	} else {
		stop("The type is not recognized")
	}
}

