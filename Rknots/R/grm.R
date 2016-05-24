# Script comments and history
# 2011
# 5:35:25 PM

# Author: Federico Comoglio @ D-BSSE, ETH Zurich
###############################################################################

findComplement <- function(set, subset) {
	if(length(set) < length(subset))
		warning("set longer than subset\n")
	complement <- set[!is.element(set, subset)]
	return(complement)
}

grm <-
function (points3D, b, ends, inputM) 
{
	M <- inputM
	n <- nrow(points3D)
	ncomp <- length(ends) + 1
	shift <- 0
	sh <- 0	
	if (b >= (n - 1) || (b %in% ends)) 
		return(list(points3D = points3D, ends = ends, M = M))
	exends <- c(0, ends, n)
	intervals <- list()
	for (k in 1 : ncomp) 
		intervals[[k]] <- (exends[k] + 1) : (exends[k + 1] - 1)
	comp <- which(exends - b >= 0)[1] - 1
	ecomponent <- exends[comp + 1]
	tempM <- M
	e <- b
	while (e < (ecomponent - 1)) {
		sh <- 0
		e <- e + 1
		range <- b : e
		usubM <- M[range, range]
		usubM[lower.tri(usubM)] <- 0
		if ((min(usubM) * max(usubM) != -1) == TRUE) {
			sh <- sh + 1
			complement <- list()
			for (k in 1 : ncomp) 
				complement[[k]] <- M[range, intervals[[k]]]
			#complement[[comp]] <- complement[[comp]][, intervals[[comp]][-which(intervals[[comp]] %in% 
			#										(b : e) == TRUE)] 
			#		- ifelse (comp > 0, exends[comp], 0)]
			complement[[comp]] <- complement[[comp]][, findComplement(intervals[[comp]], (b : e))
							- ifelse (comp > 0, exends[comp], 0)]
			
			complementsigns <- unique(unlist(complement))
			if(length(complementsigns) == 0) 
				cond <- 1
			else
				cond <- min(complementsigns) * max(complementsigns)
			if ((cond != -1) == TRUE) {
				sh <- sh + 1
				temp <- edgeIntersections(points3D, c(1 : n)[-((b + 1) : e)], b)
				nr <- rep(0, n - 1)
				nr[temp[, 1]] <- temp[, 2]
				nrtemp <- nr
				nrtemp[c(b, e)] <- 0
				nrsigns <- unique(nrtemp)
				cnd <- c(complementsigns, nrsigns)
				if ((min(cnd) * max(cnd) != 
							-1) == TRUE) {
					sh <- sh + 1
					tempM <- M
					tempM[b, ] <- nr
					tempM[, b] <- (-nr)
					tempM <- tempM[-((b + 1):e), -((b + 1):e)]
				}
				else 	break
			}
			else 	break
		}
		else	break
	}
	if (sh == 3) 
		shift <- 1
	drop.n <- shift + min(e - 1, ecomponent - 2)
	if ((b + 1) <= drop.n) {
		points3D <- points3D[-c((b + 1) : drop.n), ]
		if(comp <= (ncomp - 1)) {
			ends[comp : (ncomp - 1)] <- ends[comp:(ncomp - 1)] - (drop.n - b)
		}
		M <- tempM
	}
	return(list(points3D = points3D, ends = ends, M = M))
}

