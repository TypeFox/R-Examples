
ghg <- function(m,v,fu){
	ghg <- 0
	inf <- NULL
	sup <- NULL
	for( i in 1:m){
		for( j in 1:m){
			if(j >= i){
				ij <- (j-1)*j/2+i
				inf <- c(inf,ij)
			}else{
				ij <- (i-1)*i/2+j
				sup <- c(sup,ij)
			}
			ghg = ghg + v[m*(m+1)/2+i]*fu[ij]*v[m*(m+1)/2+j]
		}
	}
	return(list(ghg=ghg,inf=inf,sup=sup))
	
}