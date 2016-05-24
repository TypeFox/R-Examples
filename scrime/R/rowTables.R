`rowTables` <-
function(x, levels=1:3, affy=FALSE, includeNA=FALSE, 
		useNN=c("not", "only", "also"), check=TRUE){
	if(affy)
		levels <- c("AA", "AB", "BB")
	if(includeNA | check){
		useNN <- match.arg(useNN)
		if(useNN != "only")
			n.nas <- rowSums(is.na(x))
		else
			n.nas <- rowSums(x=="NN")
		if(useNN == "also")
			n.nas <- n.nas + rowSums(x=="NN", na.rm=TRUE)
	}
	if(includeNA && all(n.nas==0)){
		warning("Since x does not contain missing values,\n",
			"no column for NAs is added to the table.")
		includeNA <- FALSE
	}
	if(includeNA){			
		n.levs <- length(levels)+1
		names.levs <- c(levels, "NA")
	}
	else{
		n.levs <- length(levels)
		names.levs <- levels
	}
	tab <- matrix(0, nrow(x), n.levs, dimnames=list(rownames(x), names.levs))
	for(i in 1:length(levels))
		tab[,i] <- rowSums(x==levels[i], na.rm=TRUE)
	if(includeNA)
		tab[,n.levs] <- n.nas
	if(check){
		tmp <- rowSums(tab)
		if(!includeNA)
			tmp <- tmp + n.nas
		if(any(tmp != ncol(x)))
			warning("At least one of the rows of x seems to contain values\n",
				"differing from the specified ones.")
	}
	tab
}

