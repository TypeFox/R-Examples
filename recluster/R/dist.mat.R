dist.mat <- function(com,pair) {
	ncom <- nrow(com)
	distmat <- matrix(nrow=ncom,ncol=ncom,0,dimnames=list(rownames(com),rownames(com)))
	st <- c(0,cumsum(seq(ncom-1,2)))+1
	end <- cumsum(seq(ncom-1,1))
	for (i in 1:(ncom-1)) distmat[i,(ncom:(seq(1,ncom)[i]))]=c(pair[end[i]:st[i]],0)
	distmat <- as.dist(t(distmat))
	return(distmat)
}
