#' GEBCO plot. Plots equidepth lines.
#' 
#' Plots lines of equal depths using a database from GEBCO.
#' 
#' 
#' @param depth A vector of the depths which we want equidept lines plotted.
#' @param col The colour of the lines, if the col vector is shorter than the
#' depth vector it is repeated.  Default is all lines black.
#' @param lty Linetype, if the lty vector is shorter than the depth vector it
#' is repeated.  Default is all lines have linetype 1.
#' @param lwd Linewidth, if the lwd vector is shorter than the depth vector it
#' is repeated. Default is all lines have linewidth 1.
#' @param depthlab A boolean variable determening whether labels should be
#' printed on equidepth lines, default is false.
#' @param depthlabcex The size of depthlabels.
#' @return None
#' @section Side Effects: Plots equidepth lines on current plot.
#' @seealso \code{\link{geoplot}}, \code{\link{geolines}}.
#' @examples
#' 
#'    geoplot()   # Set up plot.
#'    
#'    gbplot(c(100,500,1000),depthlab=T,depthlabcex=0.2)
#'    # Plot depthlines for 100,500,1000,1500 m, showing the
#'    # depth on the line.   
#' 
#' @export gbplot
gbplot <-
function(depth, col, lty, lwd, depthlab, depthlabcex)
{
	if(missing(depthlabcex))
		depthlabcex <- 0.7
	if(missing(lwd))
		lwd <- rep(1, length(depth))
	if(missing(lty))
		lty <- rep(1, length(depth))
	if(missing(col))
		col <- rep(1, length(depth))
	if(length(col) < length(depth))
		col[(length(col) + 1):length(depth)] <- col[length(col)]
	if(length(lwd) < length(depth))
		lwd[(length(lwd) + 1):length(depth)] <- lwd[length(lwd)]
	if(length(lty) < length(depth))
		lty[(length(lty) + 1):length(depth)] <- lty[length(lty)]
	for(i in 1:length(depth)) {
		dypi <- depth[i]
		if(dypi %% 100 != 0 || dypi == 300 || dypi == 700) {
			print(paste(dypi, "m does not exist in GEBCO data"))
			return(invisible())
		}
		if(dypi <= 1000 || dypi == 1200 || dypi == 1500 || dypi == 2000
			)
			txt <- paste("geolines(gbdypi.", dypi, 
				",col=col[i],lwd=lwd[i],lty=lty[i])", sep = "")
		else {
			j <- match(dypi, names(geo::gbdypi))
			txt <- paste("geolines(gbdypi[[", j, 
				"]],col=col[i],lwd=lwd[i],lty=lty[i])", sep = 
				"")
		}
		eval(parse(text = txt))
		if(!missing(depthlab)) {
			k <- !is.na(match(geo::depthloc$z, dypi))
			if(any(k))
				geotext(geo::depthloc[k,  ], z = geo::depthloc[k, "z"],
					cex = depthlabcex)
		}
	}
	return(invisible())
}

