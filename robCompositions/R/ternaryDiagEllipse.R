#' Adds tolerance ellipses to a ternary diagram.
#' 
#' Low-level plot function which add tolerance ellipses to a high-level plot of
#' a ternary diagram.
#' 
#' 
#' @param x Three-part composition. Object of class \dQuote{matrix} or
#' \dQuote{data.frame}.
#' @param tolerance Determines the amount of observations with Mahalanobis
#' distance larger than the drawn ellipse, scaled to one.
#' @param locscatt Method for estimating the mean and covariance.
#' @param \dots Additional arguments passed trough.
#' @return no values are returned.
#' @author Peter Filzmoser, Matthias Templ
#' @seealso \code{\link{ternaryDiag}}
#' @keywords aplot
#' @export
#' @examples
#' 
#' data(coffee)
#' x <- coffee[,2:4]
#' ternaryDiag(x, grid=FALSE)
#' ternaryDiagEllipse(x)
#' ## or directly:
#' ternaryDiag(x, grid=FALSE, line="ellipse")
#' 
ternaryDiagEllipse <- function(x, tolerance=c(0.9,0.95,0.975), locscatt="MCD", ...){
	z <- isomLR(x)
	if(locscatt=="MCD"){
		cv <- robustbase::covMcd(z)
		mu <- cv$center
		cm <- cv$cov
	} else {
		mu <- colMeans(z)
		cm <- cov(z)
	}
	dat1 <- drawMahal(z, mu, cm, plot=FALSE, whichlines=tolerance) 
	for(i in 1:length(tolerance)){
		e <- isomLRinv(cbind(dat1$mdX[,i], dat1$mdY[,i]))
		xp1 <- e[, 2] + e[, 3]/2
		yp1 <- e[, 3] * sqrt(3)/2	  
		lines(xp1, yp1, xlim = c(0, 1), ylim = c(0, 0.9), #frame.plot = FALSE, 
				xaxt = "n", yaxt = "n", xlab = "", ylab = "", ...)
	}
}
