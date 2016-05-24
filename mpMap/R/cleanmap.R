#' Clean map for use in QTL mappin
#' 
#' Given an mpcross object, this function will remove markers which are clustered too tightly together to be useful for QTL mapping. Markers will be removed from the map, from the data, and from the estimated recombination fractions, effectively subsetting the cross down to a more spaced out grid of markers. 
#' @export
#' @param mpcross Object of class \code{mpcross}
#' @param mindist Minimum distance between markers in cM
#' @return An mpcross object is returned which has markers removed which are within the minimum distance specified. These markers do not provide additional information for QTL mapping and increase the computational burden. 
#' @seealso \code{\link[mpMap]{mpcross}}

cleanmap <- function(mpcross, mindist=1)
{
	rm <- vector()
	newmpc <- mpcross
	nmiss <- clean(mpcross)$missing
	
	offsets = c(0, cumsum(unlist(lapply(mpcross$map, length))))
	
	toRemove = c()
	for (chromosome in 1:length(mpcross$map))
	{
		map = mpcross$map[[chromosome]]
		index=1
		toRemoveThisChromosome = c()
		while(index <= length(map))
		{
			j = index + 1
			while(j <= length(map) && abs(map[j] -map[index]) < mindist) j = j + 1
			if(j != index + 1)	toRemoveThisChromosome = c(toRemoveThisChromosome, (index+1):(j-1))
			index = j
		}
		if(length(toRemoveThisChromosome) > 0) newmpc$map[[chromosome]] <- newmpc$map[[chromosome]][-1*toRemoveThisChromosome]
		toRemove = c(toRemove, toRemoveThisChromosome + offsets[chromosome])
	}
	newmpc$finals <- mpcross$finals[,-toRemove]
	newmpc$founders <- mpcross$founders[,-toRemove]
	return(newmpc)
}
