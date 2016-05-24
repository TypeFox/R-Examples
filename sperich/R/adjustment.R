adjustment <-
function(species.richness, noninterpolatedgrid, clusterlist){
	#set up sub-function
	adjust.one.cluster <- function(species.richness, noninterpolatedgrid, cluster){
		if (length(cluster)==0){
			return(species.richness)
		}

		maxc <- ifelse(max(noninterpolatedgrid[cluster]) < 0, 0, max(noninterpolatedgrid[cluster]))

		for (pixel in cluster){
			pixel.value <- ifelse(noninterpolatedgrid[pixel] < 0, 0, noninterpolatedgrid[pixel])
			quotient <- ifelse(is.nan(pixel.value/maxc),0,pixel.value/maxc)
			quotient <- ifelse(is.infinite(quotient),0,quotient)
			species.richness[pixel] <- (species.richness[pixel]-noninterpolatedgrid[pixel]) * (1 - quotient) + noninterpolatedgrid[pixel]
		}
		return(species.richness)
	}	

	#adjust species.richness for each cluster in clusterlist
	for (i in 1:length(clusterlist)){
		cluster <- clusterlist[[i]]
		species.richness <- adjust.one.cluster(species.richness, 
				noninterpolatedgrid, cluster)
	}

	return(species.richness)
}
