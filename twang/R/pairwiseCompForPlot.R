pairwiseCompForPlot <- function(mnps, stop.method = 1, returnMax = TRUE, ...){
	nFt <- mnps$nFits
	unwHldDif <- wHldDif <- matrix(NA, nrow = length(mnps$psList[[1]]$desc$unw$bal.tab$results$tx.mn), ncol = choose(nFt,2))
	cnt <- 1
	for(i in 1:(nFt-1)) for(j in (i+1):nFt){
	unwHldDif[,cnt] <- abs(mnps$psList[[i]]$desc$unw$bal.tab$results$tx.mn - mnps$psList[[j]]$desc$unw$bal.tab$results$tx.mn)/mnps$psList[[1]]$desc$unw$bal.tab$results$ct.sd
	cnt <- cnt+1
	}
	
	unwHldDif <- abs(unwHldDif)

	cnt <- 1
	for(i in 1:(nFt-1)) for(j in (i+1):nFt){
	wHldDif[,cnt] <- abs(mnps$psList[[i]]$desc[[stop.method + 1]]$bal.tab$results$tx.mn - mnps$psList[[j]]$desc[[stop.method + 1]]$bal.tab$results$tx.mn)/mnps$psList[[1]]$desc[[stop.method + 1]]$bal.tab$results$ct.sd
	cnt <- cnt+1
	}	
	
	wHldDif <- abs(wHldDif)
	
	if(!returnMax) return(list(weightedPairwiseDiff = wHldDif, unweightedPairwiseDiff = unwHldDif))
	else return(data.frame(unweightedPairwiseDiff = apply(unwHldDif, 1, max), weightedPairwiseDiff = apply(wHldDif, 1, max), row.names=row.names(mnps$psList[[1]]$desc$unw$bal.tab$results)))
	
}