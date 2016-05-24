# Author: Robert J. Hijmans, r.hijmans@gmail.com
# 2009
# Version 0.1
# Licence GPL3
   
.chao <- function(x) {
	spp <- as.vector(table(x))
	Sobs <- length(spp)
	singles <- length(spp[spp==1])
	doubles <- length(spp[spp==2])
	chao <- Sobs + singles^2/(2*doubles)
	return(chao)
}

