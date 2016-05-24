threshOpt <- function(distobj, sppVector, threshold = 0.01){
	distobj <- as.matrix(distobj)
	#individual (diag) values excluded
	diag(distobj) <- NA
	#indices of singletons
	singletons <- rmSingletons(sppVector, exclude=FALSE)
	#gets vector of column rows of all non-singletons
	ZZ <- rmSingletons(sppVector, exclude=TRUE)
	#run the sensitivity loop - singleton species (MU) excluded from iteration - only uses ZZ for query
	#False positive -  no matches within x% of query --- as the "NA" value
	#Positive - only 1 sp. within x% of query --- as the "TRUE" value;
	#False negative - more than 1 spp. within x% of query --- as the "FALSE" value
	#remember for 0, we need to say ==0 rather than <0
	OUT <- NULL
	for(i in 1:dim(distobj)[1]) {
		inThresh <- sppVector[which(distobj[,i] < threshold)]
			if(length(inThresh) == 0 && i %in% singletons) OUT[i] <- "True neg" else {
				if(length(inThresh) == 0 && !i %in% singletons) OUT[i] <- "False pos" else{
					if(length(unique(inThresh)) == 1 && unique(inThresh) == sppVector[i]) OUT[i] <- "True pos" else{
						if(length(unique(inThresh)) == 1 && !unique(inThresh) == sppVector[i]) OUT[i] <- "False neg" else OUT[i] <- "False neg"}
	}}}
	OUT <- factor(OUT, levels=c("True neg", "True pos", "False neg", "False pos"))
	tab <- table(OUT)
	tab <- c(threshold, tab, sum(tab[3], tab[4]))
	names(tab)[c(1,6)] <- c("Threshold", "Cumulative error")
	tab
}
