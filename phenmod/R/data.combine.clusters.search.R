data.combine.clusters.search <- 
		function(stations, stations.net, 
			shuffle=TRUE){
	clusters <- list()

	if (shuffle){
		stations.mod <- stations[sample(1:length(stations), length(stations), replace=FALSE)]
	} else {
		stations.mod <- stations
	}

	for (station in stations.mod){
		station.nr <- which(station==stations)
		if (length(station.nr)==0){
			next
		}
		near.stations <- stations.net[[station.nr]]
		if (is.na(near.stations[1])){
			clusters[[length(clusters)+1]] <- station
			stations[[station.nr]] <- NULL
			stations.net[[station.nr]] <- NULL
		} else {
			clusters[[length(clusters)+1]] <- c(station, near.stations)
			for (station.to.delete in c(station, near.stations)){
				station.to.delete.nr <- which(station.to.delete==stations)
				stations[[station.to.delete.nr]] <- NULL
				stations.net[[station.to.delete.nr]] <- NULL
				if (length(stations.net)>0){
					for (nr in 1:length(stations.net)){
						next.near.stations <- stations.net[[nr]]
						if (length(which(next.near.stations == station.to.delete)) > 0){
							next.near.stations <- next.near.stations[which(next.near.stations != station.to.delete)]
							if (length(next.near.stations) > 0){
								stations.net[[nr]] <- next.near.stations
							} else {
								stations.net[[nr]] <- NA
							}
						}
					}
				}
			}
		}
	}
	
	return(clusters)
}