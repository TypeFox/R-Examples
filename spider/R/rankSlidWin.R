rankSlidWin <- 
function(slidWin, criteria = "mean_distance", num = 10){
	revRank <- function(xx) (length(xx)+1) - rank(xx, ties.method="min")
	list2df <- function(listObj){
		len <- sapply(listObj, length)
		mlen <- min(len)
		for(i in 1:length(listObj)) listObj[[i]] <- listObj[[i]][1:mlen]
		as.data.frame(listObj)
	}
	distCriteria <- distLabel <- treeCriteria <- treeLabel <- treeEX <- NULL
	if(!slidWin$distMeasures && !slidWin$treeMeasures) stop("Object of class`slidWin' must be created using slideAnalyses")
	if(!slidWin$distMeasures && slidWin$treeMeasures) criteria <- "monophyly"
	if(slidWin$distMeasures){
		distCriteria <- c("mean_distance", "zero_noncon", "zero_distances", "diag_nuc")
		distLabel <- c("dist_mean_out", "noncon_out", "zero_out", "nd_out")
	}
	if(slidWin$treeMeasures){ 
		treeCriteria <- c("monophyly", "clade_comparison", "clade_comp_shallow")
		treeLabel <- c("win_mono_out", "comp_out", "comp_depth_out")
		treeEX <- c("pos_tr_out")
	}
	measures <- c("position", distCriteria, treeCriteria)
	#Remove objects not of interest
	excluded <- match(c("dat_zero_out", "boxplot_out", "distMeasures", "treeMeasures", treeEX), names(slidWin))
	dFrame <- list2df(slidWin[-excluded])
	#Reorder and rename dataframe columns
	dfOrder <- match(c("pos_out", distLabel, treeLabel), names(dFrame))
	dFrame <- dFrame[ , dfOrder]
	names(dFrame) <- measures
	#Order rows
	high <- match(c("monophyly", "clade_comparison", "clade_comp_shallow", "diag_nuc", "mean_distance"), measures)
	highVal <- as.data.frame(lapply(dFrame, revRank))
	dfVals <- highVal
	if(slidWin$distMeasures){
		low <- match(c("zero_noncon", "zero_distances"), measures)
		lowVal <- as.data.frame(lapply(dFrame, function(x) rank(x, ties.method="min")))
		dfVals[ , low] <- lowVal[ , low]
	}
	if("all" %in% criteria) rowOrd <- order(apply(dfVals[ , -1], MARGIN=1, FUN=sum))
		else{
			ordNum <- which(measures %in% criteria)
			if(length(criteria) > 1) rowOrd <- order(apply(dfVals[ , ordNum], MARGIN=1, FUN=sum))
				else rowOrd <- order(dfVals[ , ordNum])
				}
	#Return top 10
	head(dFrame[ rowOrd , ], n = as.integer(num))
}