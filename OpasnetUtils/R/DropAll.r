# DROPALL #####################################
## dropall pudottaa data.framesta pois kaikki faktorien sellaiset levelit, joita ei kayteta.
## parametrit: x = data.frame   

dropall <- function(x){
	isFac <- NULL
	for (i in 1:dim(x)[2]){
		isFac[i] = is.factor(x[ , i])
	}
	for (i in 1:length(isFac)){
		x[, i] <- x[, i][ , drop = TRUE]
	}
	return(x)
}
