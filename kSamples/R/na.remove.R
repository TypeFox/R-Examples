na.remove <- function(x){
#
# This function removes NAs from a list and counts the total 
# number of NAs in na.total.
# Returned is a list with the cleaned list x.new and with 
# the count na.total of NAs.
#
	na.status <- lapply(x,is.na) # changed sapply to lapply
	k <- length(x)
	x.new <- list()
	na.total <- 0
	for( i in 1:k ){
		x.new[[i]] <- x[[i]][!na.status[[i]]]
		na.total <- na.total + sum(na.status[[i]])
	}
	list(x.new=x.new,na.total=na.total)
}
