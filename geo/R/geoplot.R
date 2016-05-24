#' Plots lat and lon coordinates using Mercator or Lambert transformation.
#' 
#' Creates a plot of given data, as the "plot" function but allows input data
#' to have latitude and longitude coordinates.Coordinates are latitude and
#' longitude in degrees (default) or x, y in m (or other units).
#' 
#' The program performs the Mercator(default) or Lambert transformation of the
#' data in degrees and plots it. The plot is scaled such that 1cm in vertical
#' and horizontal directions correspond to the same distance. The program is
#' used to initialize drawings to be used by programs, i.e. geolines,
#' geopolygon, geopoints, geotext, geosymbol, geogrid, geocontour.fill,
#' geoimage and geocontour. The inputs to the program are decimal numbers
#' negative for western longitudes e.g.lat = 65.75 lon = -25.8. Whether the
#' program interprets data as lat, lon or x, y depends on the parameter
#' projection. If projection is "none" data is interpreted as x, y else lat,
#' lon. When the data is interpreted as x, y lists are assumed to have the
#' components \$x and \$y, else \$lat and \$lon.Geoplot is often used with type
#' = "n" so the datapoints used to set up the drawing are not seen. Calling
#' geoplot with plotit = FALSE sets up the drawing without putting anything on
#' the screen ( or page).
#' 
#' 
#' @param lat,lon Latitude and longitude of data (or x and y coordinates),
#' negative for southern latitudes and western longitudes. May be supplied as
#' two vectors or as a list or dataframe lat (or x) including vectors lat\$lat
#' and lat\$lon (x\$x and x\$y if projection = none). If xlim and ylim are
#' given the arguments lat and lon are not required.
#' @param type Options are the same as in the plot function i.e "l" for lines,
#' "n" for not plotting the data, "p" for points etc. Default value is "p"
#' @param pch Type of symbol drawn in points. Default is "*". Any other
#' character or symbol can be used, seet the help on points.
#' @param xlim x limits of drawing (longitudinal direction). If not given the
#' program finds the limits as the range of the data times the parameter r
#' (1.05 default).  xlim can be a list (dataframe) with components lat and lon
#' (x and y if projection = "none").
#' @param ylim y limits of drawing (latitudinal direction). If not given and
#' xlim is not a list the program finds the ulimits as the range of the data
#' times the parameter r (1.05).
#' @param b0 Base latitude for the Mercator or Lambert transform.  Default
#' value is 65 (typical for Iceland).  If Mercator transform is used values
#' close to the mean of the data are recommended except the data extend very
#' far north, (> 75 th degree) then b0 should be close to the northern limits
#' of the data).
#' @param r Size of area which is plotted. r = 1.0 means exactly the range of
#' the datapoints, but r = 1.5 means that the range is 1.5 times the range of
#' data. Default value is 1.05. Does not matter if xlim and ylim are given.
#' @param country Country that is plotted.  Options are: \verb{ 1. island 1300
#' points (iceland) 2. bisland 20000 points (iceland fine) 3. greenland 62000
#' points 4. faeroes 2100 points 5. eyjar 2200 points(islands around iceland)
#' 6. none no country plotted } Which map is used depends on the size of the
#' area plotted. If small part of the coast is seen bisland is used, else
#' Iceland. Isands can be added later by geolines(eyjar). Default is set by a
#' variable with the name COUNTRY.DEFAULT.  COUNTRY.DEFAULT <- "island" makes
#' iceland the default country.  The program geoworld can be used to add all
#' coastlines that intersect the map and also to make new countries (see
#' geoworld).
#' @param xlab X-label. Default value is " "
#' @param ylab Y-label. Default value is " "
#' @param option Can be either "cut" or "nocut". If "nocut" the plot always
#' fills the plotted area but if "cut" the plot does not fill it in the
#' direction where the range of data is minimum. It has to be kept in mind that
#' the program always keeps the same scale vertically and horizontally. Default
#' value is "cut". Not effective when contourplots are plotted.
#' @param grid If grid is TRUE meridians and paralells are plotted, else not.
#' Default value is TRUE.
#' @param new If new is FALSE the plot is added to the current plot otherwise a
#' new plot is made. Default value is FALSE. Similar to the Splus command
#' "par(new = TRUE)".
#' @param cont A parameter to reserve space for legend for contours besides the
#' plot. Default value is "FALSE".  Rarely used at the legends are usually put
#' somewhere in the graph.
#' @param cex Relative size of character and symbols (see the help on the
#' parameter cex).  The size of plotted characeters is cex time the parameter
#' csi that can be seen by par()\$csi.  In earlier versions of geoplot the
#' parameter csi was set but csi is a parameter that can not be set in R.
#' @param col The color number used to for the plot. Default value is 1 that is
#' usually black.
#' @param lcont Limits of area preserved for lables and contour plot as ration
#' of plotting area. Default value is c(0.13, 0.21). That means that the labels
#' take 13\% of the plotting area but the figure 79\%. (Labels to left, figure
#' to left.) 0.08 to 0.1 is a reasonable difference. Only used when cont = TRUE
#' which is seldom used as described erlier.
#' @param plotit If FALSE plot is only initialized but not plotted. If used
#' other programs are used to fill the plot (geolines, geocontour, geopolygon
#' etc). Most often used in multiple plots.
#' 
#' Used in connection with geocontour.fill to fewer files but the plot command
#' is given again with new = TRUE when geocontour.fill is called. Plot = FALSE
#' does not work if axeslabels = FALSE. Something seems to have to be on the
#' graph for proper setup.
#' @param reitur Reitur means statistical square in Icelandic.  If true the
#' division of the axes is idendical to the distribution of the ocean around
#' Iceland in statistical squares. Means dlat = 0.5;dlon = 1.
#' @param smareitur If true the division of the axes is idendical to the
#' distribution of the ocean around Iceland in subsquares.dlat = 0.25;dlon =
#' 0.5
#' @param reittext If true the number of each square is written in the center
#' of the square. Only meaningful for Icelandic waters.
#' @param cexrt Relative size of reittext text in statistical squares
#' (cexrt*csi). Default value is 0.7.
#' @param axratio Parameter usually not changed by the user.
#' @param lwd Line width for plot (grid and axes). Default value is the value
#' set when the program was called.(usually 1). Higher values correspond to
#' wider lines.
#' @param lwd1 Line width for plot country. Default value is the value set when
#' the program was called (usually 1). Higher values correspond to wider lines.
#' @param locator Some kind of a simple zoom command. If locator is TRUE the
#' user points with the mouse on the limits of the plot he wants to make. Needs
#' a plot made by geoplot on the screen. Same as if zoom is used.
#' @param axlabels If FALSE no numbers are plotted on the axes. Default value
#' is TRUE.
#' @param projection Projection used to make the plot. Options are "Mercator",
#' "none" and "Lambert". Default value is "Mercator". If projection = "none"
#' data is assumed to be x, y.
#' @param b1 Second latitude to define Lambert projection.
#' @param dlat Defines the grid, to make a grid on the lat axis, 1 is a number
#' on axis and a line at every deg. Not usualy set by user.
#' @param dlon Same as dlat, but for lon.
#' @param jitter useful if many datapoints have the same coordinates, points
#' are jittered randomly to make common values look bigger. jitter = 0.016 is
#' often about right but you may want to have jitter smaller or bigger varying
#' on plot.
#' @param zoom If TRUE the mouse is used to point to two points on the current
#' map with those two new points becoming the new corner points on a new map.
#' . Same as the parameter locator in older geoplot versions.
#' @param csi Size of character.  This parameter can not be set in R but for
#' compatibility with old Splus scripts the parameter cex is readjusted by cex
#' = cex*csi/0.12.  Use of this parameter is not recommended.  Default value is
#' NULL i.e not used.
#' @param csirt Size of reittext. Default value is 0.1. Only for compatibility
#' with older Splus programs but cexrt should be used as the parameter csi can
#' not be set in R:
#' @param xaxdist Distance from plot to the labels on the xaxis (dist or r
#' argument to geoaxis.  Default value is 0.2 but higher value mean that
#' axlabels is further away from the plot.  Further flexibility with axes can
#' be reached by calling geoplot with axlabels = FALSE and geoaxis aferwards.
#' @param yaxdist Distance from plot to the labels on the yaxis (dist or r
#' argument to geoaxis.  Default value is 0.3 but higher value mean that
#' axlabels is further away from the plot.
#' @return No values are returned. The graphical setup is stored in a global
#' list called geopar. That list is accessed by other program that use the same
#' setup.
#' @section Side Effects: There should be no side effects. The program changes
#' a number of graphical parameters but the old parameters are restored before.
#' @seealso \code{\link{geolines}}, \code{\link{geopolygon}},
#' \code{\link{geotext}}, \code{\link{geosymbols}}, \code{\link{geogrid}},
#' \code{\link{geopar}}, \code{\link{geocontour.fill}},
#' \code{\link{geolocator}}, \code{\link{geocontour}},
#' \code{\link{reitaplott}}, \code{\link{geodefine}}, \code{\link{Proj}}.
#' @examples
#' 
#' \dontrun{Examples shown here below also include calls to the other functions
#' in the geopackage.  Further explanations of these functions can be
#' found in the appropriate help files.
#' 
#' Contour plot of haddock catch in icelandic waters based on logbooks.
#' A color scheme
#' where color 0 is white, 1 black and 2-150 gradually changing from
#' white to black is used.  The data is for 6 years and is stored in a
#' list hadcatch with 6 components.  Circles showing haddock catch in
#' the Icelandic groundfish survey is added on top of the plot with the
#' function geosymbols but utbrteg is a dataframe with information on all
#' catch in the Icelandic groundfish survey (all names in Icelandic ysa
#' means haddock and ar year).  The function bwps is a call to the
#' postscript function with the indicated color scheme.  Designed for
#' Splus and has to be changed for R as the color schemes there are quite
#' different.  This applies to all the examples below.
#' 
#' lev <- c(0.5, 1, 2, 4, 6)
#' col <- c(0, 30, 50, 70, 90, 150)
#' txt <- c(1993, 1995, 2000, 2002, 2005, 2006)
#' par(mfrow = c(3, 2));par(mex = 0.01)
#' bwps(file = "hadcatchutbr.ps", height = 6.8, width = 6.5, horizontal = FALSE)
#' for(i in 1:6) {
#'  SMB.std.background(grid = FALSE, axlabels = FALSE)
#'  geocontour.fill(hadcatch[[i]], levels = lev, col = col,
#'     white = TRUE, working.space = 2e6)
#'  gbplot(200)
#'  geotext(67.3, -27.6, txt[i], csi = 0.16, adj = 0)
#'  geopolygon(island, col = 0);geolines(island)
#'  tmp <- utbrteg[utbrteg$ar == txt[i], ] # select the year
#'  geosymbols(tmp, z = tmp$ysa.kg, circles = 0.2,
#'     sqrt = TRUE, lwd = 1)# amount of haddock
#'  geopoints(tmp, pch = 16, csi = 0.05)
#' }
#' dev.off()
#' 
#' # plot x, y data
#' geoplot(x$x, x$y, projection = "none", type = "n")
#' 
#' geoplot(x, projection = "none", type = "n")
#' # does the same thing.
#' 
#' # The packages maps and mapdata need to be installed
#' # worldHires is a very detailed database of coastlines from the
#' # package mapdata.  Could be problematic if used with fill = TRUE)
#' # Allowed.size is the maximum allowed size of polygons.
#' library(map) # world coastlines and programs
#' library(mapdata) # more detailed coastlines
#' geoplot(xlim = c(20, 70), ylim = c(15, 34))
#' geoworld(database = "worldHires", fill = TRUE, col = 30, allowed.size = 30000)
#' 
#' geoplot(xlim = c(20, 70), ylim = c(15, 34), dlat = 10, dlon = 10)
#' geoworld(database = "world", fill = TRUE, col = 30) #
#' 
#' geoplot(xlim = c(-10, 70), ylim = c(71, 81), b0 = 80,
#'   dlat = 2, dlon = 10) # 0 must be high here else
#' geoworld(database = "world", fill = TRUE, col = 30) #the plot fails.
#' 
#' # Lambert projection,
#' geoplot(xlim = c(-10, 70), ylim = c(71, 81),
#'   dlat = 2, dlon = 10, projection = "Lambert")
#' geoworld(database = "world", fill = TRUE, col = 30)
#' 
#' # Lambert projection, get the axis closer with the mgp command
#' par(mgp = c(2, 0, 0))
#' geoplot(xlim = c(-10, 70), ylim = c(71, 81),
#'   dlat = 2, dlon = 10, projection = "Lambert", cex = 1.1)
#' geoworld(database = "world", fill = TRUE, col = 30)
#' 
#' # Example with capelin data.  (lodna meanns capelin).  lodna.2 is a
#' # data.frame with components lat,  lon and z.
#' geoplot(lodna.2, type = "l")
#' geopoints(lodna.2)
#' geosymbols(lodna.2, z = lodna.2$z, colplot = TRUE,
#'   parbars = 0.05, levels = vor.levels, label.location = labloc)
#' 
#' 
#' geoplot(lodna.2, type = "l")
#'  geosymbols(lodna.2, z = lodna.2$z, perbars = 0.1)
#' 
#' geoplot(lodna.2, type = "l")
#'  geosymbols(lodna.2, z = log(1+lodna.2$z), perbars = 0.1)
#' 
#' 
#' limits <- list(lat = c(63, 68), lon = c(-30, -10))
#' geoplot(xlim = limits, type = "n", grid = FALSE, axlabels = TRUE, plot = FALSE)
#' tmp <- geoexpand(lodna.2.grd) # expand the grid
#' # has defined ther area vor.area and data outside it are set  to NA.
#' i <- geoinside(tmp, vor.area.new, option = 0)
#' zgr <- z.lodna.2;zgr[-i] <- NA
#' geocontour.fill(lodna.2.grd, z = zgr, white = TRUE,
#'   label.location = labloc, levels = vor.levels)
#' geoplot(xlim = limits, type = "n", grid = FALSE, axlabels = TRUE, new = TRUE)
#' #geolines(lodna.2, lwd = 1)
#' gbplot(c(200, 500)) # Depth contours.
#' 
#' # make a plot of number within a square, calculate the total number of
#' #cod (torskur) within a square (reitur) and put the text number of
#' #square and total number of cod in the center of the square (number of
#' #cod below number of square).  Apply.shrink is similar to tapply
#' #returning the data in different form and is included with the geo library
#' 
#' geoplot(island, r = 1.2, type = "n", reitur = TRUE)
#' x <- apply.shrink(data$torskur.stk, data$reitur, sum,
#'   names = c("reitur", "torskur.stk"))
#' x1 <- r2d(x$reitur)
#' geotext(x1, z = paste(x$reitur, round(x$torskur.stk, 1), sep = "\n"))
#' 
#' # Plot filled circles.  The color scheme used is the same as described
#' # color 0 white, 1 black and 2 - 155 white-black see bwps
#' # in geosymbols the argument color means size (in inches)
#' # when the fill.circles = TRUE.  The data used  AfliBySquareMonthYear have
#' # the columns year, month , square and catch.
#' # the function r2d changes square (reitur in Icelandic) to position
#' # Text is put in the middle of the circles where catch exceeds 1000
#' # tonnes.
#' 
#' my.colors = c(0.004, 0.04, 0.1, 0.15, 0.20, .25, 100)
#' lev <- c(0.2, 2, 7.5, 10, 20, 50)
#' 
#' yy <- c(1932:1939)
#' tmp4 <- AfliBySquareMonthYear
#' tmp4$catch <- tmp4$catch/1000
#' for (ar in yy) {
#'   bwps(file = paste(ar, ".ps", sep = ""))
#'   par(omi = c(0, 0, 0, 2))
#'   par(mfrow = c(4, 3))
#'   par(mex = 0.01)
#' 
#'   for(man in 1:12){
#'     tmp1 <- tmp4[tmp4$year == ar & tmp4$month == man, ]
#'     SMB.std.background(axlabels = FALSE, country = "none", plotit = FALSE)
#'     if(nrow(tmp1) > 0) {
#'       tmp2 <- apply.shrink(tmp1$catch, tmp1$square, sum,
#'         names = c("square", "catch"))
#'       tmp2 <- tmp2[order(-tmp2$catch), ]
#' 
#'       tmp3 <- data.frame(r2d(tmp2$square))
#'       geosymbols(tmp3, z = tmp2$catch, fill.circles = TRUE, col = 60,
#'         levels = lev, colors = my.colors, bordercol = 0, border = TRUE)
#'       geopolygon(island, col = 30)
#'       geolines(eyjar, lwd = 3, col = 30)
#'       j <- tmp2$catch > 1
#'       if(any(j)) {
#'         tmp3 <- tmp3[j, ]
#'         tmp2 <- tmp2[j, ]
#'         geotext(lat = tmp3$lat, lon = tmp3$lon, z = tmp2$catch,
#'           angle = 45, csi = 0.1)
#'       }
#'       geotext(lat = c(65.2), lon = (-18), z = paste(month.abb[man],
#'         round(sum(tmp2$catch, na.rm = TRUE)), sep = "\n"), csi = 0.2)
#'     }
#'     else {geotext(lat = c(65.2), lon = (-18),
#'             z = paste(month.abb[man], "0", sep = "\n"), csi = 0.2)}
#'   }
#'   geotext(63, -10, ar, adj = 1, csi = 0.18)
#'   dev.off()
#' }
#' }
#' @export geoplot
geoplot <-
function(lat = NULL, lon = 0, type = "p", pch = "*", xlim = c(0, 0), 
  ylim = c(0, 0), b0 = 65, r = 1.05, country = "default", xlab = " ", 
  ylab = " ", option = "cut", grid = TRUE, new = FALSE, cont = FALSE,
  cex = 0.9,col = 1, lcont= c(0.13, 0.21), plotit = TRUE, 
  reitur = FALSE, smareitur = FALSE, reittext = FALSE, cexrt = 0.7, 
  csirt=NULL, axratio = 1, lwd = 0, lwd1 = 0, locator = FALSE, 
  axlabels = TRUE, projection = "Mercator", b1 = b0, dlat = 0, dlon = 0, 
  jitter = 0, zoom, csi = NULL, xaxdist = 0.2, yaxdist = 0.3)
{
  geopar <- getOption("geopar")
  if(!is.null(csirt)) cexrt <- cexrt*csirt/0.12
  if(!is.null(csi)) cex <- cex*csi/0.12
  if(!plotit) axlabels <- FALSE  # not plot axes if ther is no plot.  
  if(!missing(zoom)) {
    xlim <- geolocator(n = 2)
  }
	oldpar.1 <- par(no.readonly = TRUE)
 	# first version of old parameters
	command <- sys.call()
	if((oldpar.1$fig[2] - oldpar.1$fig[1]) <= 0.6 || (oldpar.1$fig[4] -
		oldpar.1$fig[3]) <= 0.6)
		multfig <- TRUE
	else multfig <- FALSE
	if(projection == "none") {
                if(is.list(xlim) &&  any(!is.na(match(c("x","y"),names(xlim))))){
                  ylim <- xlim$y
                  xlim <- xlim$x
		}
	}
	else {
                if(is.list(xlim) &&  any(!is.na(match(c("lat","lon"),names(xlim))))){
                  ylim <- xlim$lat
                  xlim <- xlim$lon
		}
	}
	if(is.null(lat) && xlim[2] == xlim[1] && ylim[2] == ylim[1] && !locator
		) {
		#std plot
		if(!multfig) par(fig = geo::geopar.std$fig)
		if(!multfig)
			par(plt = geo::geopar.std$plt)
		xlim <- geo::geopar.std$xlim
		ylim <- geo::geopar.std$ylim
		if(!multfig)
			par(mex = geo::geopar.std$mex)
	}
	if(is.null(lat)) {
		lat <- c(65, 66)
		lon <- c(-28, -27)
		type <- "n"
	}
	oldpar <- selectedpar()
	if(locator & missing(zoom)) {
		limits <- geolocator(n = 2)
		if(geopar$projection == "none") {
			xlim <- limits$x
			ylim <- limits$y
		}
		else {
			xlim <- limits$lon
			ylim <- limits$lat
		}
	}
	xlim <- sort(xlim)
	ylim <- c(ylim)
	if(projection == "none") {
		if(length(country) == 1)
			if(country == "default")
				country <- "none"
	}
	else {
		if(length(country) == 1)
                  if(country == "default")
                    eval(parse(text = paste("country <- ",geo::COUNTRY.DEFAULT)))
	}
	init(lat, lon = lon, type = type, pch = pch, xlim = xlim, ylim = ylim,
		b0 = b0, r = r, xlab = xlab, ylab = ylab,
		option = option, grid = grid, new = new, cont = cont, cex = cex,
		col = col, lcont = lcont, plotit = plotit, reitur = reitur,
		smareitur = smareitur, reittext = reittext, axratio = axratio,
		lwd = lwd, axlabels = axlabels, oldpar = oldpar, projection = 
		projection, b1 = b1, dlat = dlat, dlon = dlon, command = 
		command, jitter = jitter, xaxdist = xaxdist, yaxdist = yaxdist)
        oldpar.1 <- Elimcomp(oldpar.1)
	par(oldpar.1)
#        par(new = TRUE)
	if(reittext)
		plot_reitnr(cexrt, lwd = lwd)
	# number of squares
	if(length(country) > 1 && plotit) geolines(country, col = col, lwd = 
			lwd1)
	# plot country
	return(invisible())
}

