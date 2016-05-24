`flexmonoreg.poolpvalues` <-
function(m, x, lambda=0, alpha.location=1, alpha.adjacency=0.5, B=100){

pvalue.temp1 <- function(t, d){
	return( sum(d>t)/length(d) )
}
pvalue.temp2 <- function(d, t){
	return( sapply(t, pvalue.temp1, d) )
}

flexmonoreg.temppoolstat <- function(k, m, x, lambda=0, alpha.location=1, alpha.adjacency=0.5){
	return( apply(m, 1, flexmonoreg.stat, sample(x), lambda, alpha.location, alpha.adjacency) )
}

fscore <- apply(m, 1, flexmonoreg.stat, x, lambda, alpha.location, alpha.adjacency)

fscore.perm <- sapply(c(1:B), flexmonoreg.temppoolstat, m, x, lambda, alpha.location, alpha.adjacency)

return( rowMeans(apply(abs(fscore.perm), 2, pvalue.temp2, abs(fscore))) )
}