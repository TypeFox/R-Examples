`flexmonoreg.pvalue` <-
function(y, x, lambda=0, alpha.location=1, alpha.adjacency=0.5, B=100){

flexmonoreg.tempstat <- function(k, y, x, lambda=0, alpha.location=1, alpha.adjacency=0.5){
	return( flexmonoreg.stat(y, sample(x), lambda, alpha.location, alpha.adjacency) )
}

fscore <- flexmonoreg.stat(y, x, lambda, alpha.location, alpha.adjacency)

fscore.perm <- sapply(c(1:B), flexmonoreg.tempstat, y, x, lambda, alpha.location, alpha.adjacency)

return( sum(fscore.perm>fscore)/B )
}