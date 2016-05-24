create_DF_rank <- function(data, ord.by, group, median.row=FALSE, rev.ord=FALSE, align='top'){
	# data=dStats; align=vertical.align; group=grouping
	DF <- data

	DF$rank <- rank(DF[,ord.by], ties.method = "first")# create a new rank column
	if(rev.ord) DF$rank <- max(DF$rank)+1 - DF$rank
	DF <- DF[order(DF$rank),]
	
	m.rank <- (max(DF$rank)+1)/2
	if(median.row) DF$median <- (DF$rank==m.rank) else DF$median <- FALSE

	group <- group[!group=='median']
	if(length(group)==1) group <- rep(group, ceiling(sum(1-DF$median)/group))

      iGroups2 <- cumsum(group) 
	if(median.row){
	  fixGrouping=FALSE
	  is.region.median <- nrow(DF)%%2
	  warning0 <- (sum(group) < (nrow(DF)-1))
	  if(warning0) stop("Grouping does not account for all data frame rows")

	  warning1 <- !(sum(group) == nrow(DF) - is.region.median)
	  warning2 <- !(sum(group[iGroups2 < m.rank]) == floor(m.rank) - is.region.median)
	    if(warning1) warning("Grouping specification does not match dataframe row count.", call. = FALSE, immediate. = TRUE)
	    if(warning2) warning("Grouping does not adequately allow for median.", call. = FALSE, immediate. = TRUE)
	    if(warning1 | warning2) fixGrouping=TRUE
	  
	  if(fixGrouping){
	    warning("Reminder -- Do not specify a median row in grouping arguement. Attempting grouping auto-alteration", call. = FALSE, immediate. = TRUE)
	    w.tophalf <- sum(iGroups2 < m.rank)
	    tmpGrouping <- NULL
	    try(tmpGrouping <- c(group[1:w.tophalf], 
			  floor(m.rank-1/2) - iGroups2[w.tophalf],
			  iGroups2[w.tophalf+1] - floor(m.rank-1/2),
			  group[-(1:(w.tophalf+1))]))
	    if(!is.null(tmpGrouping)){
		tmpGrouping <- tmpGrouping[!tmpGrouping==0]
		group <- tmpGrouping 
	    	iGroups2 <- cumsum(group) 
	    }
	    if(is.null(tmpGrouping)) stop("Auto-alteration of groupings failed", call. = FALSE)
	  }
	  if(is.region.median) iGroups2 <- c(iGroups2[iGroups2<m.rank], max(iGroups2[iGroups2<m.rank]) + 1, iGroups2[iGroups2>m.rank]+1)
	}
	
	DF$pGrp <- as.numeric(cut(DF$rank, c(0,iGroups2), labels=1:length(iGroups2))) 
	    # create a new perceptual group column based on rank
		
			
	pGrpStats <- aggregate(DF$rank, list(DF$pGrp, DF$median), length)
	names(pGrpStats) <- c('pGrp','median','length')	
	pGrpStats$addOrd <- 0
	if(align=='center') pGrpStats$addOrd <- (1-pGrpStats$median)*(max(pGrpStats$length)-pGrpStats$length)/2
	
	DF <- merge(DF, pGrpStats[,c('pGrp','addOrd')])
	
	DF$pGrpRank <- sapply(1:nrow(DF), function(i) sum(DF$rank[i]>subset(DF, DF$pGrp==DF$pGrp[i])$rank)+1)
	DF$pGrpOrd <- DF$pGrpRank + DF$addOrd
	DF$color <-DF$pGrpRank
		
	DF
}


create_DF_cat <- function(data, grp.by, cat){
	DF <- data 

	tGroups <- unique(DF[,grp.by])
	DF$pGrp <- match(DF[,grp.by],tGroups)									
	tCats <- unique(DF[,cat])
	DF$pGrpOrd <- match(DF[,cat], tCats)	
	DF$color <-DF$pGrpOrd

	DF								
}



alterForMedian <- function(DF, a){
	if(a$median.row){ 
	    if(!any(DF$median)) {
		tmpCols <- names(DF)[-(1:a$ncols)]
		tmpData <- apply(DF[,tmpCols],2,median)

		DFmedian <- transform(DF[1,], pGrpOrd=1, pGrp=a$m.pGrp, median=TRUE, rank='')
		for(k in 1:length(tmpCols)) DFmedian[,tmpCols[k]] <- tmpData[k]

		DF <- rbind(DF, DFmedian)
	      }

	    DF$color[DF$median] <- length(a$colors)
	  }

	DF
}
