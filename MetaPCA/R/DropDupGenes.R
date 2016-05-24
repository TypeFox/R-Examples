#When multiple probesets share the same gene symbols, select only the best probeset in terms of IQR
DropDupGenes <- function(dat, isParallel=FALSE, nCores=NULL, na.rm=TRUE) {
	.dropped <- if(na.rm) which(is.na(rownames(dat)) | rownames(dat)=='') else vector()
	.g <- factor(rownames(dat))
	.dupIDs <- levels(.g)[table(.g)>1]
	
	if(length(.dupIDs)>0) {
		if(isParallel) {
			.workers <- NULL
			if(!is.null(nCores))
				options(cores=nCores)
			if(.Platform$OS.type == "unix") {
				requireAll("doMC")
				registerDoMC()
			} else { #windows
				requireAll("doSMP")
				.workers <- startWorkers() #default is 3
				registerDoSMP(.workers)
			}
			.dropped <- c(.dropped, foreach(i=iter(.dupIDs), .combine=c) %dopar% {
						which.g <- which(.g==i)
						which.g[-which.max(apply(dat[which.g,],1,IQR,na.rm=T))]
					})
			if(!is.null(.workers))
				stopWorkers(.workers)
		} else {
			for(i in .dupIDs) {
				which.g <- which(.g==i)
				.dropped <- c(.dropped, which.g[-which.max(apply(dat[which.g,],1,IQR,na.rm=T))])
			}
		}
		return(dat[-.dropped,])
	} 
	if(length(.dropped)>0)
		return(dat[-.dropped,])
	else
		return(dat)
}
