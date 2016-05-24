uvector <- function(v){

	if(sum(is.na(v)) > 0) return(v)

	if(is.vector(v)){
		d <- sqrt(sum(v^2))
		if(d == 0){return(c(0,0,0))}
		return(unlist(v/d))
	}
	if(is.matrix(v)){
		d <- sqrt(apply(v^2, 1, sum))
		for(i in 1:length(d)){
			if(d[i] == 0) next
			v[i, ] <- v[i, ]/d[i]
		}
		return(v)
	}
}