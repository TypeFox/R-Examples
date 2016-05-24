quade <- function(obj){
	Qs <- lapply(obj,function(st) st$rank*as.vector(st$possibleTreat %*% st$withinRank))
	Qtab <- do.call(expand.grid,Qs)
	Q <- rowSums(Qtab)
	Qobs <- sum(sapply(obj,function(st) st$rank*as.vector(st$obsTreat %*% st$withinRank)))
	Qprobs <- apply(do.call(expand.grid,lapply(obj, function(st) st$prob)),1,prod)
	pval <- sum(Qprobs[which(Q>=Qobs)])
	perms <- data.frame(Q=Q,prob=Qprobs)
	out <- list(Qobs=Qobs,Q=Q,permutations=perms,pval=pval,sets=obj)
	attr(out,"pairs") <- attr(obj,"pairs")
	class(out) <- "quade"
	return(out)
}


print.quade <- function(x,...){
	if(attr(x,"pairs")){
		cat(paste("Observed signed rank statistic (Qobs):",x$Qobs,"\n",sep=" "))
	} 
	else{
		cat(paste("Observed Quade's statistic (Qobs):",x$Qobs,"\n",sep=" "))
	}
	cat(paste("Pr(>= Qobs):",round(x$pval,5),"\n",sep=" "))
}