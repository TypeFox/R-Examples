#' Function that fits a model to a variogram.
#' 
#' Function that fits a model to a variogram.  The fitting occurs either
#' automatically or interactively.  The function is called after the function
#' variogram which calculates the variogram to which the model is fitted.
#' Currently only spherical model is supported.  Later other models will be
#' added.
#' 
#' 
#' @param vagram List with the calculated variogram. Components of the list
#' are: \$dist mean distance of the interval.  <s-example> \$vario calculated
#' value of the variogram.  \$number number of datapoints in the interval.
#' </s-example> In nearly all cases vagram will be the output from the program
#' variogram.
#' @param model Type of model.  Default is spherical.  It is currently the only
#' model supported.
#' @param option Method to use in automatic fitting.  Allowed values are 1,2,3
#' and 4.  Default value is 2.  For further information see below.
#' @param interactivt If T the fitting is done interactively by plotting the
#' variogram on the screen and asking the user to select sill, range and nugget
#' by the locator function.
#' @param sill Sill of the variogram, or: Limit of the variogram tending to
#' infinity lag distances (wikipedia).
#' @section Value: <s-example> A list with the following components.  \$nugget
#' : Estimated nugget effect \$sill : Estimated sill \$range : Estimated range
#' \$dist : mean distance of the interval.  \$vario : calculated value of the
#' variogram.  \$number : number of datapoints in the interval. </s-example>
#' @seealso \code{\link{variogram}}, \code{\link{pointkriging}}.
#' @export variofit
variofit <-
function(vagram, model = 1, option = 2, interactivt = F, sill = 0)
{
	if(model == 1 && !interactivt && option < 5) {
		if(sill == 0 && length(vagram$variance) > 0)
			vagram$variance <- sill
		vgr <- fitspher.aut.1(vagram, option, sill)
		if(vgr$error == 1)
			return()
		return(list(rang1 = vgr$rang1, sill = vgr$sill, nugget = vgr$
			nugget, dist = vagram$dist, vario = vagram$vario, 
			number = vagram$number))
	}
	if(interactivt) {
		k <- 1
		ld <- floor(length(vagram$dist)/1.2)
		sill <- rang1 <- nugget <- rep(0, 10)
		ans <- "y"
		col <- 10
		plvar(vagram, fit = F)
		#		cat(" What type of model : , Gaussian, spherical ")
		while(ans == "y" || ans == "Y") {
			cat(" Give sill, range and nugget  in this order : \n")
			x <- locator(n = 3)
			sill[k] <- x$y[1]
			rang1[k] <- x$x[2]
			nugget[k] <- x$y[3]
			txt0 <- paste(" k = ", as.character(k))
			txt1 <- paste(" sill =", as.character(round(sill[k],
				digits = 2)))
			txt2 <- paste("range = ", as.character(round(rang1[
				k], digits = 2)))
			txt3 <- paste("nugget = ", as.character(round(nugget[
				k], digits = 2)))
			print(txt1)
			print(txt2)
			print(txt3)
			lines(vagram$dist, spherical(rang1[k], sill[k], nugget[
				k], vagram$dist), col = col)
			text(vagram$dist[ld], sill[k], as.character(k), col = 
				col)
			col <- col + 10
			cat(" \n Try again y/n  : ")
			ans <- readline()
			k <- k + 1
			if(k == 10)
				break()
		}
		nm <- 0
		cat("Give the number of the best model, default the last one:")
		nm <- scan(n = 1)
		if(length(nm) == 0)
			nm <- k - 1
		# default.  
		rang1 <- rang1[nm]
		sill <- sill[nm]
		nugget <- nugget[nm]
	}
	return(list(rang1 = rang1, sill = sill, nugget = nugget, dist = vagram$
		dist, vario = vagram$vario, number = vagram$number))
}

