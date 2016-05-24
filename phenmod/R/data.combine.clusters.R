data.combine.clusters <- function(dataset, stations.net, 
					shuffle=TRUE, tries=100, 
					silent=FALSE, out2File=FALSE){
	# get all stations
	stations <- as.list(unique(dataset$STAT_ID))
	
	# output message preparation
	msg <- ""
	
	# prepare cluster list
	clusters <- list()

	# iterate over all tries
	for (try in 1:tries){
		clusters.tmp <- data.combine.clusters.search(
					stations, stations.net, shuffle)
		if ((try==1)||(length(clusters.tmp) < length(clusters))){
			clusters <- clusters.tmp
		}
		rm(clusters.tmp)
		
		## output message
		if (!silent){
			if (out2File){
				cat("\n")
			} else {
				cat(rep("\b",nchar(msg)),sep="")
			}
			msg <- paste(try, " of ", tries, " tries done!",sep="")
			cat(msg,sep="")
			flush.console()
		}
	}
	if (!silent){
		cat("\n",sep="")
	}

	return(clusters)
}