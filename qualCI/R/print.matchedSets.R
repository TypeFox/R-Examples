print.matchedSets <- function(x,...){
	nsets <- length(x)
	sets <- names(x)
	if(attr(x,"unitNames")){
		tab <- t(sapply(x,function(st) c(paste(names(which(st$obsTreat==1)),collapse=", "),paste(names(which(st$obsTreat==0)),collapse=", "), st$rank)))
		tab <- cbind(sets,tab)
		colnames(tab) <- c("Set", "Treated","Control","Rank")
	}
	else{
		tab <- t(sapply(x,function(st) c(sum(st$obsTreat==1),sum(st$obsTreat==0), st$rank)))
		tab <- cbind(sets,tab)
		colnames(tab) <- c("Set", "Num. Treated","Num. Control","Rank")
	}
	if(attr(x,"pairs")){
		cat(paste("Design contains",nsets,"matched pairs:\n\n",sep=" "))
	}
	else{
		cat(paste("Design contains",nsets,"matched sets:\n\n",sep=" "))	
	}
	print(as.data.frame(tab),right=FALSE,quote=FALSE, row.names=FALSE)
}