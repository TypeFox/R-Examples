idQuan <- function(h,q){
	if(is.na(h)) return(NA)
	
	pos <- as.numeric(which(sort(c(h,q))==h)[1])
	return(pos)
}