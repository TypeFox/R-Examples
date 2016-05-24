## Generic function for extracting balance tables from ps and other objects
bal.table.hidden <- function(x, pairwise = TRUE, digits = 3){
	if(class(x) != "mnps"){
   bal.tab <- lapply(x$desc, function(x){return(round(x$bal.tab$results, digits))})
   return(bal.tab)
   }
   else {
   	nFits <- x$nFits
   	balTabList <- vector(mode = "list", length = nFits)
   	if(x$estimand == "ATT")
   	cat(paste("Note that `tx` refers to the category specified as the treatATT, ", x$treatATT, ".\n\n", sep = ""))
   	for(i in 1:nFits) balTabList[[i]] <- bal.table(x$psList[[i]], digits = digits)
   	if(x$estimand == "ATT") names(balTabList) <- x$levExceptTreatATT
   	if(x$estimand == "ATE") names(balTabList) <- x$treatLev
   	if(pairwise){
   		allDiffs <- NULL
   		for(i in 1:length(balTabList[[1]])){
   			holdDiffs <- rep(0, nrow(balTabList[[1]][[1]]))
   			for(j in 2:nFits){
   				for(k in 1:(j-1)){
   					holdDiffs <- apply(cbind(holdDiffs, abs(balTabList[[j]][[i]]$tx.mn - balTabList[[k]][[i]]$tx.mn)), 1, max, na.rm=TRUE)
   				}
   			}
   			allDiffs <- cbind(allDiffs, holdDiffs/balTabList[[1]][[1]]$ct.sd)
   		}
   		row.names(allDiffs) <- row.names(balTabList[[1]][[1]])
   		dimnames(allDiffs)[[2]] <- names(x$psList[[1]]$desc)
#   		return(list(pairwiseDiff = allDiffs, balanceTable = balTabList))
   		return(balTabList)
   		
   	}
   	else return(balTabList)
   }
   	
   	
}


