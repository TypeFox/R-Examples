#' Plots a variogram estimated by variofit.
#' 
#' Plot a variogram and a spherical model fit to it.
#' 
#' 
#' @param vagram List from the program variogram with the following components.
#' <s-example> \$range: range of data.  \$nugget: nugget effect.  \$sill: sill.
#' \$vario: z values of the variogram.  \$xh: mean distance for each class.
#' 
#' </s-example>
#' @param n The plot covers n times the range.  Default value is all the data.
#' @param fit Also plot lines and fitted values (\code{TRUE}, the default, 
#' or plot variogram only as points (\code{FALSE}). 
#' @param type Plot type, given to \code{plot}, defaults to points.
#' @return No values returned.
#' @seealso \code{\link{pointkriging}}, \code{\link{grid}}.
#' @export plvar
plvar <-
function(vagram, n = 4, fit = T, type = "p")
{
	if(fit) {
		vagr1 <- vagram$dist[vagram$dist < vagram$rang1 * n]
		# n x range
		zva <- vagram$vario[1:length(vagr1)]
		plot(vagr1, zva, xlim = c(0, max(vagr1)), ylim = c(0, max(
			zva) * 1.05), xlab = "Distance", ylab = "Variogram",
			title = " ", type = type)
		lines(c(0, vagr1[1]/2, vagr1), spherical(vagram$rang1, vagram$
			sill, vagram$nugget, c(0, vagr1[1]/2, vagr1)))
		tloc <- c(max(vagr1) * 0.9, max(zva) * 1.04)
		tloc <- matrix(tloc, 2, 2, byrow = T)
		tloc[2, 2] <- max(zva) * 1.01
		tmp <- vagram$nugget/vagram$sill
		tmp <- round(tmp, digits = 2)
		tmp <- as.character(tmp)
		tmp <- substring(tmp, 1, 4)
		txt <- c(paste("nugget/sill=", tmp), paste("range = ", 
			as.character(round(vagram$rang1, digits = 2))))
		print(txt)
		text(tloc, txt)
	}
	else {
		plot(c(0, vagram$dist), c(0, vagram$vario), ylim = c(0, max(
			vagram$vario)), xlab = "Distance", ylab = "Variogram",
			title = " ")
	}
}

