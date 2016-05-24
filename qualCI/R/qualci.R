qualCI <- function(obj){
	# if matched pairs
	if(attr(obj,"pairs")){
		lowest.rank.set <- which.min(sapply(obj$sets,function(x) x$rank))
		lowest.rank.set.name <- names(lowest.rank.set)
		# order by treated unit first
		lowest.rank.pair <- names(sort(obj$sets[[lowest.rank.set]]$obsTreat, decreasing=TRUE))
	}
	else{
		candidatePairsBySet <- lapply(obj$sets, function(st) getCandidatePairs(st))
		candidatePairsBySet <- candidatePairsBySet[sapply(candidatePairsBySet,length)!=0]
		candidatePairsDF <- do.call(rbind ,lapply(1:length(candidatePairsBySet), function(s) data.frame("Pair"=candidatePairsBySet[[s]], "Set"=names(candidatePairsBySet)[s])))
		cat("All treatment / control pairs of adjacent rank:\n")
		print(candidatePairsDF)
		cat("\n")
		hardpair.idx <- 0
		while(hardpair.idx < 1 | hardpair.idx > nrow(candidatePairsDF)){
			hardpair.idx <- readline("Please enter the number of the pair that was hardest to rank and press enter: ")
			hardpair.idx <- ifelse(grepl("\\D", hardpair.idx),0,as.integer(hardpair.idx))	
			if(is.na(hardpair.idx)){stop("You have stopped the execution of qualci function.",call.=FALSE)}
		}
		lowest.rank.set.name <- candidatePairsDF[hardpair.idx,"Set"]
		lowest.rank.pair <- strsplit(as.character(candidatePairsDF[hardpair.idx,"Pair"])," / ")[[1]]
	}
	# ci level
	alpha <- within.statistic(obj$sets)
	ci.level <- 1-alpha
	out <- list(qualci=lowest.rank.pair, set=lowest.rank.set.name, ci.level=ci.level)
	attr(out,"pairs") <- attr(obj,"pairs")
	class(out) <- "qualCI"
	return(out)
}


print.qualCI <- function(x,...){
	if(attr(x,"pairs")){
		cat(paste("Lower bound of qualitative one-sided confidence interval based on sign test can be described as the difference between ",x$qualci[1], " and ", x$qualci[2], " in set ", x$set, "\n", sep=""))
		cat(paste("Confidence level: ", round(x$ci.level*100,2), "%", sep=""))
	} else{
		cat(paste("Lower bound of qualitative one-sided confidence interval based on stratified rank sum test can be described as the difference between ",x$qualci[1], " and ", x$qualci[2], " in set ", x$set, "\n", sep=""))
		cat(paste("Confidence level: ", round(x$ci.level*100,2), "%", sep=""))
	}
}


getCandidatePairs <- function(st){
	rank.order <- order(st$withinRank, decreasing=T)
	treat.ordered <- st$obsTreat[rank.order]
	breaks <- (treat.ordered[-1L] != treat.ordered[-length(treat.ordered)]) & (treat.ordered[-length(treat.ordered)]==1)
	ind.breaks <- as.numeric(which(breaks))
	out <- sapply(ind.breaks, function(i) paste(names(treat.ordered[i]),"/",names(treat.ordered[i+1]),sep=" "))
	out
}

