#' Adds a plot to an existing plot initialized by geoplot.
#' 
#' Adds a plot, of any kind, to the current plot at a location specified by the
#' user.
#' 
#' The \pkg{geo} functions rely on the parameters in \code{options("geopar")}
#' for plotting. Every time \code{geoplot} is called a new set of geopar is
#' made and the current erased. If you want to change the parameters in the
#' background plot after you call geosubplot you have to save the current
#' geopar before calling geosubplot. When reassigning geopar, the command
#' \code{options(geopar=list(...))} must be used.
#' 
#' When using \code{geosubplot} it is important to save these parameters before
#' the subplot is plotted if one wants to make changes to the large plot
#' afterwards.
#' 
#' @param fun a call to the function to be subplotted.  Be sure to call geoplot
#' with new = T.
#' @param pos Position of the plot.  Should include \$lat and \$lon, if
#' length(\$lat) = 1 that point will determine the middle of the subplot, if
#' length(\$lat)=2 the points will determine oppesite corners of the subplots
#' position.
#' @param size The width and the height of the plot if length(pos\$lat)=1.
#' Default is 2,2.
#' @param fill If true the background of the subplot area will be filled.
#' Default is F.
#' @param fillcol fill colour, the colour of the background of the subplot
#' area, default is 0 (usually white)
#' @param \dots The program accepts any parameter that the subplot program will
#' accept.
#' @return none
#' @section Side Effects: A plot added to the current plot.
#' @seealso \code{\link{subplot}}, \code{\link{geoplot}}, \code{\link{geopar}}.
#' @examples
#' 
#' \dontrun{
#'       #####################################################
#'       # Example 1                                         #
#'       #####################################################
#' 
#'       geoplot()
#'       pos<- list(lat = 63.5,lon=-11.75)
#'       geosubplot(geoplot(faeroes,pch=" ",country =faeroes,new=T)
#'                  ,pos,fill=T)
#'       # Plots the Faeroes on a plot with Iceland. Be sure to
#'       # use new = T if geoplot is called again.
#' 
#'       #####################################################
#'       # Example 2                                         #
#'       #####################################################
#' 
#'       geoplot()
#'       large.geopar <- geopar             # Parameters saved.
#'       pos <- list(lat=c(63,64),lon=c(-27,-24))
#'       geosubplot(geoplot(island, new=T,grid=F,type="l"),pos)
#' 
#'       geotext(65,-18,"subplot")           # Text on subplot.
#'       small.geopar <- geopar # Parameters for subplot saved.
#' 
#'       # geopar <- large.geopar       # Make big plot active.
#'       # Unless you are working directly with the .Data dir of the
#'       # geolibrary this assignment will not work, must use:
#'       assign("geopar",large.geopar,where=0)
#' 
#'       geotext(65,-18,"Big plot")         # Text on big plot.
#' 
#'       # Another subplot.
#' 
#'       pos <- list(lat=c(63,64),lon=c(-17,-14))
#'       geosubplot(geoplot(island,new=T,grid=F,type="l"),pos,fill=T)
#'       # Another subplot drawn.
#'       geotext(65,-18,"subplot # 2")
#' 
#'       small.geopar.2 <- geopar # parameters for subplot # 2 saved.
#'       # geopar <- large.geopar # Big plot made active again
#'       # Same as above, instead use:
#'       assign("geopar",large.geopar,where=0)
#' 
#'       # See also similar example in geopar.
#' }
#' @export geosubplot
geosubplot <-
function(fun, pos, size = c(2, 2), fill, fillcol, ...)
{
	geopar <- getOption("geopar")
	if(length(pos$lat) == 1) {
		# Calculate new limits.
		plt.size <- par()$pin
		rlon <- (diff(geopar$origin$lon) * size[1])/plt.size[1]
		rlat <- (diff(geopar$origin$lat) * size[2])/plt.size[2]
		pos <- data.frame(lat = pos$lat + c(-0.5, 0.5) * rlat, lon = 
			pos$lon + c(-0.5, 0.5) * rlon)
	}
	if(!missing(fill)) {
		if(!missing(fillcol))
			geopolygon(pos, col = fillcol)
		else geopolygon(pos, col = 0)
	}
	pos <- Proj(pos)
	oldpar <- selectedpar()
	par(geopar$gpar)
	on.exit(par(oldpar))
	pr <- subplot(fun, pos, ...)
	return(invisible())
}

