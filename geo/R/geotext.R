#' Plots text on a drawing defined by geoplot.
#' 
#' The function plots a text in each of the points defined by lat,lon.  (or
#' x,y) The value plotted at each point is defined by the numeric vector or
#' text vector z. Simular to text but allows lat and lon and a choice of digits
#' used.
#' 
#' 
#' @param lat,lon Latitude and longitude of data ( or x and y coordinates),
#' negative for southern latitudes and western longitudes.  May be supplied as
#' two vectors or as a list lat (or x) including vectors lat\$lat and lat\$lon
#' (x\$x and x\$y if projection = none).
#' @param z Vector with values that will be plotted at datapoints. Has to be of
#' the same length as lat. If z is of mode character it is written directly on
#' the screen.
#' @param cex Relative size of characters (see the help on the parameter cex).
#' The size of plotted characeters is cex time the parameter csi that can be
#' seen by par()\$csi.  In earlier versions of geoplot the parameter csi was
#' set but csi is a parameter that can not be set in R.
#' @param adj Location of the text relative to the point, 0 means text right of
#' point, 0.5 text centered at point and 1 text left of point.  Default value
#' is 0.5.
#' @param col Color number used. Default value is 1.
#' @param digits Number of digits used.
#' @param pretext Text put in front of all the text.  Default value is nothing.
#' @param lwd Linewidth used.
#' @param aftertext Text put after all the text.  Default value is nothing.
#' @param outside If outside is F no text is plotted outside the range defined
#' by xlim,ylim.  Else it is done.  Default value is F.
#' @param angle angle of text in degrees, default is 0.
#' @param jitter see jitter in geoplot.
#' @param csi Size of character.  This parameter can not be set in R but for
#' compatibility with old Splus scripts the parameter cex is readjusted by cex
#' = cex*csi/0.12.  Use of this parameter is not recommended.  Default value is
#' NULL i.e not used.
#' @return No values returned.
#' @seealso \code{\link{geoplot}}, \code{\link{geopolygon}},
#' \code{\link{geolines}}, \code{\link{geosymbols}}, \code{\link{geogrid}},
#' \code{\link{geopar}}, \code{\link{geocontour.fill}},
#' \code{\link{geolocator}}, \code{\link{geocontour}}.
#' @examples
#' 
#' \dontrun{       geotext(deg,z=z)    # plot text at points deg$lat,deg$lon
#'        
#'        geotext(deg$lat,deg$lon,z,csi=0.06) # Same, size of text 0.06".
#'        
#'        geotext(x$x,x$y,x$z,aftertext="km",pretext="distance")
#'        # If geopar$projection="none"
#'        
#'        geotext(x,z=x$z,aftertext=" km",pretext="distance",angle =90) 
#'        # Same text written vertically.
#' 
#' 
#'        ###############################################################
#'        #  Example                                                    #
#'        ###############################################################
#' 
#'        lon <- rnorm(10,-27,1.3)
#'        lat <- rnorm(10,65,0.6)
#'        # Make a normal dist. random set of 10 points.
#' 
#'        geoplot(lat=lat,lon=lon,grid=F,xlim=c(-22,-30),ylim=c(63,67))
#'        # Plot the random data points.
#'        
#'        geopolygon(island,col=115,exterior=T)
#'        geolines(island)
#'        # Color Iceland. Use litir(number) to see colour scheme. 
#'        # Sharpen lines around Iceland.
#' 
#'        num <- 1:10
#'        lab <- paste("Nr.",num,sep="")
#'        # Make string vector with "Nr.1".."Nr.10" for geotext.
#' 
#'        geopoints(lat,lon,pch="*",col=5)
#'        # Redraw the data in a new color the * mark at points.
#' 
#'        geotext(lon=lon,lat=lat,z=lab,col=155)
#'        # With geotext we put one element from lab at each data point.
#'        title(main="10 Random Data Point")
#'        # Add title
#' }
#' @export geotext
geotext <-
function(lat, lon = 0, z, cex = 0.7, adj = 0.5, col = 1, digits = 0, pretext
	 = "", lwd = 0, aftertext = "", outside = F, angle = 0, jitter = NULL,csi=NULL)
{
    geopar <- getOption("geopar")
    if(!is.null(csi)) cex <- cex*csi/0.12 # For compatibility
	if(length(lon) == 1 && length(lat) > 1) {
		if(geopar$projection == "none") {
			lon <- lat$y
			lat <- lat$x
		}
		else {
			lon <- lat$lon
			lat <- lat$lat
		}
	}
	if(geopar$projection != "none") {
		# degrees and minutes
		if(mean(lat, na.rm = T) > 1000) {
			lat <- geoconvert(lat)
			lon <-  - geoconvert(lon)
		}
	}
	if(!is.null(jitter)) {
		lat <- lat + runif(length(lat), -1, 1) * jitter
		lon <- lon + runif(length(lon), -1, 1) * jitter * 2
	}
        
	oldpar <- selectedpar()
	par(geopar$gpar)
        if(outside)
		par(xpd = T)
	else par(xpd = F)
	if(lwd != 0)
		par(lwd = lwd)
	on.exit(par(oldpar))
	par(cex = cex)
	par(adj = adj)
	xx <- Proj(lat, lon, geopar$scale, geopar$b0, geopar$b1, geopar$l1,
		geopar$projection)
	ind <- c(1:length(xx$x))
	if(is.character(z)) {
		if(pretext == "")
			txt <- z
		else txt <- paste(pretext, z, sep = "")
		if(aftertext != "")
			txt <- paste(txt, aftertext, sep = "")
	}
	else {
		if(pretext == "")
			txt <- format(round(z, digits = digits))
		else txt <- paste(pretext, format(round(z, digits = digits)),
				sep = "")
		if(aftertext != "")
			txt <- paste(txt, aftertext, sep = "")
	}
	if(!outside) {
		ind <- c(1:length(xx$x))
		ind <- ind[(!is.na(xx$x)) & (xx$x < geopar$limx[1] | xx$x >
			geopar$limx[2] | xx$y < geopar$limy[1] | xx$y > geopar$
			limy[2])]
		xx$x[ind] <- NA
		xx$y[ind] <- NA
	}
	if(length(angle) == length(xx$x) || length(col) == length(xx$x)) {
		if(length(angle) < length(xx$x))
			angle <- rep(angle[1], length(xx$x))
		if(length(col) < length(xx$x))
			col <- rep(col[1], length(xx$x))
		for(i in 1:length(xx$x)) {
			text(xx$x[i], xx$y[i], txt, col = col[i], srt = angle[
				i])
		}
	}
	else text(xx$x, xx$y, txt, col = col, srt = angle)
	return(invisible())
}

