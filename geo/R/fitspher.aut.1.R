#' Fits a super-smoother to a variogram (?)
#' 
#' Fits a super-smoother to a variogram (?).
#' 
#' 
#' @param vagram Variogram
#' @param option What to do ?
#' @param sill Sill for the variogram (?)
#' @return Returns a list with components of the returned variogram:
#' \item{nugget}{Nugget} \item{dist}{Distance} \item{range}{Range}
#' \item{sill}{Sill} \item{error}{Does the function call return an error? Used
#' in \code{variofit}}
#' @note Needs further elaboration.
#' @seealso Called by \code{\link{variofit}}.
#' @keywords smooth
#' @export fitspher.aut.1
fitspher.aut.1 <-
function(vagram, option, sill)
{
	i1 <- c(2:(length(vagram$vario) - 1))
	# index vector
	zvar1 <- supsmu(vagram$dist, vagram$vario)
	# Do the fitting, find range, sill and nugget effect.
	# Max value of supersmoother.  		
	if(option == 1) ind <- i1[zvar1$y[i1] == max(zvar1$y[i1])]
	# First maximum of supersmoother
	if(option == 2) ind <- i1[(zvar1$y[i1] > zvar1$y[i1 + 1]) & (zvar1$
			y[i1] > zvar1$y[i1 - 1])][1]
	# Second maximum of supersmoother.   
	if(option == 3) ind <- i1[(zvar1$y[i1] > zvar1$y[i1 + 1]) & (zvar1$
			y[i1] > zvar1$y[i1 - 1])][2]
	# Sill given find range
	if(option == 4) {
		if(sill == 0)
			sill <- vagram$variance
		ind <- i1[zvar1$y[i1] > sill][1]
	}
	if(is.na(ind) && (option != 1)) {
		print(" condition not satisfied.  Max value of")
		print(" supersmoother used ")
		ind <- i1[zvar1$y[i1] == max(zvar1$y[i1])]
		option == 1
	}
	if(length(ind) == 0 && (option == 1)) {
		print(" cannot fit variogram")
		error <- 1
		return(error)
	}
	rang1 <- zvar1$x[ind]
	# range
	if(option != 4) sill <- zvar1$y[ind]
	# sill
	xvar2 <- (1.5 * vagram$dist[1:ind])/rang1 - (0.5 * vagram$dist[1:ind]^
		3)/rang1^3
	zvar2 <- vagram$vario[1:ind] - sill * xvar2[1:ind]
	xvar2 <- 1 - xvar2
	nugget <- lsfit(xvar2, zvar2,  , F)$coef
	# 	if the nugget effect is estimated less than 0 it set to 0.05  
	if(nugget < 0) {
		nugget <- 0.05
	}
	error <- 0
	return(list(nugget = nugget, dist = dist, range = rang1, sill = sill,
		error = error))
}

