rank_na <- function (dat, na=-100000){
## Transform the original matrix in a rank matrix where the NAs have been replaced by the value na.
	n <- nrow(dat)
	m <- ncol(dat)
	for(i in 1:m){
		dat[,i][is.na(dat[,i])==FALSE]<- rank(dat[,i][is.na(dat[,i])==FALSE])
	}
	dat[which(is.na(dat)==TRUE)] <- na
	dat <- matrix(dat,n,m)
	return(dat)
}

