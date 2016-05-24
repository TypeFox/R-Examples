`flexisoreg.poolpvalues` <-
function(m, x, lambda=0, alpha.location=1, alpha.adjacency=0.5, B=100){

pvalue.temp1 <- function(t, d){
	return( sum(d>t)/length(d) )
}
pvalue.temp2 <- function(d, t){
	return( sapply(t, pvalue.temp1, d) )
}

flexisoreg.temppoolstat <- function(k, m, x, lambda=0, alpha.location=1, alpha.adjacency=0.5){
	return( apply(m, 1, flexisoreg.stat, sample(x), lambda, alpha.location, alpha.adjacency) )
}

fscore <- apply(m, 1, flexisoreg.stat, x, lambda, alpha.location, alpha.adjacency)

fscore.perm <- sapply(c(1:B), flexisoreg.temppoolstat, m, x, lambda, alpha.location, alpha.adjacency)

return( rowMeans(apply(fscore.perm, 2, pvalue.temp2, fscore)) )
}