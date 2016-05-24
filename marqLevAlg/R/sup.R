sup <- function(m){
	indice<-NULL 
	for(i in 1:m){
		for(j in i:m){
			indice <- c(indice,(j-1)*j/2+i)
		}
	}
	return(indice)
}