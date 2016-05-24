# TODO: Add comment
# 
# Author: ian30
###############################################################################


## Rekursive Funktion zur Berechnung der Kombinationen der Anzahl
## der Erkrankten in den Risikogruppen:
computeDiseaseMatrix <- function(nDis, maxVec){
	if (length(maxVec) == 1){
		out <- matrix(nDis, ncol = 1)		
	} else {
		from <- max(c(0,(nDis-sum(maxVec[-1]))))
		to <- min(c(nDis,maxVec[1]))
		if (from > to){
			out <- matrix()
		} else {
			firstCol <- max(c(0,(nDis-sum(maxVec[-1])))):min(c(nDis,maxVec[1]))		
			out <- NULL
		    for (ii in seq(along = firstCol)){
				out <- rbind(out, 
					cbind(firstCol[ii], computeDiseaseMatrix(nDis-firstCol[ii], 
					maxVec[-1])))				
			}			
		}		
	}
	return(out)
}
