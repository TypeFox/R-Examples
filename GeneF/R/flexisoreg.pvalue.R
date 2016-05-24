`flexisoreg.pvalue` <-
function(y, x, lambda=0, alpha.location=1, alpha.adjacency=0.5, B=100){

flexisoreg.tempstat <- function(k, y, x, lambda=0, alpha.location=1, alpha.adjacency=0.5){
	return( flexisoreg.stat(y, sample(x), lambda, alpha.location, alpha.adjacency) )
}

fscore <- flexisoreg.stat(y, x, lambda, alpha.location, alpha.adjacency)

fscore.perm <- sapply(c(1:B), flexisoreg.tempstat, y, x, lambda, alpha.location, alpha.adjacency)

return( sum(fscore.perm>fscore)/B )
}