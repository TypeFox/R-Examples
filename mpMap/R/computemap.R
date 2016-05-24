#' Computes map distances
#' 
#' Given an mpcross object with a map order and matrix of recombination fractions, this function will estimate map positions
#' @export
#' @param object Object of class \code{mpcross}
#' @param mapfx Map function to convert recombination fractions to cM
#' @param missfx Function to impute missing recombination fraction values
#' @return An mpcross object is returned whose map component has been estimated based on the map order and matrix of recombination fractions. Missing recombination fractions are imputed either by filling in the closest non-missing value (missfx=1) or by averaging the distance between other nearby markers (missfx=2). 
#' @seealso \code{\link[mpMap]{mpcross}}


computemap <- function(object, mapfx=c("haldane", "kosambi"), missfx=2)
{
	if (missing(mapfx)) mapfx <- "haldane"
	if (mapfx=="haldane") mf <- haldaneR2X else mf <- kosambiR2X

	if (is.null(object$rf)) object <- mpestrf(object)

	if (is.null(object$map)) {
		cat("Warning: no map order listed. Assuming all markers are on the same chromosome\n")
		object$map <- list()

		rf <- fill(fill(object$rf$theta, missfx), 1)
          
		rf[rf==.5] <- .49
		object$map[[1]] <- c(0, cumsum(mf(rf[row(rf)==(col(rf)+1)])))
		names(object$map[[1]]) <- colnames(object$finals)
	} else {
		for (i in 1:length(object$map)){
			m <- match(names(object$map[[i]]), colnames(object$finals))
		
	   		rf <- fill(fill(object$rf$theta[m,m], missfx), 1)
            
			rf[rf==.5] <- .49
			object$map[[i]] <- c(0, cumsum(mf(rf[row(rf)==(col(rf)+1)])))
			names(object$map[[i]]) <- colnames(object$finals)[m]
		}
	}
	return(object)
}


