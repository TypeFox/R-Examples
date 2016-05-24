classToInd <- function(x){
	
	count.vec <- vector()
	for (i in 1:nrow(x)){
		tmp <- rep(x[i,1], x[i,2])
		count.vec <- c(count.vec, tmp)
	}
	ind.table <- cbind(c(1:length(count.vec)), count.vec)
	colnames(ind.table) <- c("ID", "No.Obs")
	return(ind.table)
	
}