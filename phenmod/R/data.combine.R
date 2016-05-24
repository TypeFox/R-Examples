data.combine <- function(dataset, range=5000, alt.range=50, shuffle=TRUE, 
				tries=100, silent=FALSE, out2File=FALSE,
				clusters.tmp.file="tmpcluster.RData"){
	# create clusters of stations
	stations.net <- data.combine.stationNet(dataset, 
				range, alt.range, silent, out2File)
	clusters <- data.combine.clusters(dataset, 
			stations.net, shuffle, tries, silent, out2File)
	if (!is.null(clusters.tmp.file)){ 
		save(clusters, file=clusters.tmp.file) 
	}

	# combine time series
	data.combined <- data.combine.timeseries(dataset, clusters, 
				silent, out2File, minimalClusterSize=5)

	return(data.combined)
}	
