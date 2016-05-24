within.statistic <- function(obj){
	Qs <- lapply(obj,function(st) as.vector(st$possibleTreat %*% st$withinRank))
	Qtab <- do.call(expand.grid,Qs)
	Q <- rowSums(Qtab)
	Qobs <- sum(sapply(obj,function(st) as.vector(st$obsTreat %*% st$withinRank)))
	Qprobs <- apply(do.call(expand.grid,lapply(obj, function(st) st$prob)),1,prod)
	pval <- sum(Qprobs[which(Q>=Qobs)])
	perms <- data.frame(Q=Q,prob=Qprobs)
	return(pval)
}