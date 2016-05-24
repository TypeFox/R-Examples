stat_Chen2010 <- function(sam1, sam2){
	n1 <- dim(sam1)[1]
	n2 <- dim(sam2)[1]
	col.sum1 <- colSums(sam1)
	col.sum2 <- colSums(sam2)
	rm1 <- sum(col.sum1^2) - sum(diag(sam1 %*% t(sam1)))
	rm2 <- sum(col.sum2^2) - sum(diag(sam2 %*% t(sam2)))
	chen2010.stat <- rm1/(n1*(n1 - 1)) + rm2/(n2*(n2 - 1)) - 2*sum(col.sum1*col.sum2)/(n1*n2)
	chen2010.stat <- as.numeric(chen2010.stat)
	return(chen2010.stat)
}