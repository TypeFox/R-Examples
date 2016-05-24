#' Plot different kinds of symbols at the data points.
#' 
#' The function plots different kinds of symbols at the data points defined by
#' lat, lon. There are four categories of symbols: 
#'
#' \describe{
#'  \item{\strong{Default}}{Shapes whose size is proportional to z or sqrt(z).}
#'  \item{\strong{Categories}}{Shapes where certain color, shading or size
#'   represents certain range of z. Similar to contour program.}
#'  \item{\strong{Filled circles}}{Certain size represents certain range of z,
#'   specified with fill.circles = TRUE.}
#'  \item{\strong{Characters}}{Characters or character strings represent the
#'   different ranges of z, specified with characters = TRUE.}
#' }
#'
#' There are seven types of shapes: circles, squares,
#' rectangles, vbars (vertical bars), hbars(horisontal bars) , perbars
#' (perpendicular bars), parbars (parallel bars)
#' 
#' 
#' @param lat,lon latitude and longitude of data or a dataframe containing
#' latitude and longitude of data (or x and y coordinates), negative for
#' southern latitudes and western longitudes. Expected to contain \$lat and
#' \$lon if not otherwise specified in col.names.
#' @param z Matrix containing values at datapoints.
#' @param levels Values at contourlines. Default value is zero. If levels is
#' zero the program determines the contourlines from the data. If squares,
#' circles, hbars, vbars or perbars are being used levels corresponds to the
#' levels where labels are given.
#' @param reflevels i Some levels for reference?
#' @param labels.only if true only labels are plotted. Default is false.
#' @param cex Size expansion of digits.
#' @param chs Something to do with characters?
#' @param z1 Value 2 at data points. Only used in connection with rectangles.
#' @param circles Max size of circles plotted at data points. Default value
#' used if value <-0 or >100. Size is either proportional to z or sqrt(z).
#' @param squares Max size of squares plotted at data points. Default value
#' used if value <-0 or >100. Size is either proportional to z or sqrt(z).
#' @param rectangles Max size of rectangles plotted at datapoints in inches.
#' The first number gives max height and the second number max width. Values <
#' 0 or > 100 give default values.
#' @param vbars Max size of vertical bars at data points in inches. Values >100
#' give default values. Value <0 gives bars below points.
#' @param hbars Max size of horizontal bars at data points in inches. Values
#' >100 give default values. Value <0 gives bars left of points.
#' @param perbars Max size of bars perpendicular to transsect lines in inches.
#' Values >100 give default values. Value <0 gives different orientation.
#' @param parbars Same as perbars, except the bars are now parallel.
#' @param sqrt If sqrt = TRUE the size of symbol is proportional to sqrt(z)
#' else to z. Default value if FALSE.
#' @param col Color number used, default col = 0.
#' @param maxn If nonzero maxn is the base for the size of symbols, else max(z)
#' is used. Nonzero maxn is used if several plots are to be compared.
#' @param colplot if true range of z is specified by a colour and not size.
#' @param nlevels Number of contourlines. Used if the program has to determine
#' the contourlines. Default value is 10.
#' @param colors Color number for the contourlines. Runs from 0 to 155. A
#' vector one longer than the vector levels. Default is blue-green-yellow- red
#' from lowest to the highest values. Color 0 is white and 1 is black. On black
#' and white plots higher color number means darker color and also when
#' hatching. When hatching the useful range is 10 to 80. When plotting filled
#' circles of different sizes colors means sizes in inches. When levels are not
#' specified directly colors have to be found by the program because even
#' though nlevels = 5 the length of levels can be 7 due to characteristics of
#' the Splus pretty command.
#' @param n Number of vertices in each circle, default value is 25, this
#' parameter is rarely changed by the user.
#' @param maxcol Number of colors used (excluding #0). Default value is 155
#' @param only.positive Logical value. If FALSE then negative values are
#' allowed else negative values are set to zero. Default value is FALSE.
#' @param digits Number of digits used in labels. Default value is zero.
#' @param white If true the first color is white.
#' @param lwd Line width for symbols. Default value is the value when the
#' program was called.
#' @param label.location List with components \$lat and \$lon specifying
#' oppesite corners of a square where the label should be put. (or \$x, \$y)
#' Gives the lower left and upper right corner of the box where the labels are
#' put. Default value is 0 that means no labels are put on the drawing or if
#' geoplot was initialized with cont = TRUE, then labels are put on the left
#' side of the plot. label.location is best given by geolocator or directly by
#' specifying label.location = "locator".
#' @param labels Type of labels 1 or 2. One is default and is usually used
#' except in color with very many colors (more than 10-20) .
#' @param fill.circles If TRUE filled circles of different sizes are plotted.
#' Size of the circles is given directly by the parameter color or found from
#' the data. (maxcol corresponds to the size of the largest circle in this case
#' and is 0.1 by default (changing maxcol to 0.4 makes all the circles 4 times
#' larger.)
#' @param density If density is 1 (or not zero) circles are hatched instead of
#' having different color. Only available with circles. Color does in this case
#' specify the density of hatching. Higher number means denser hatching. The
#' range is from zero to maxcol (155). But the effective range is ca. 10 - 80.
#' @param angle Angle of hatching, default is 45 degrees.
#' @param rotate Rotation of hatching from one level to the next. Default value
#' is 0 but 45 or 90 can be good to better distinquish between different
#' levels.
#' @param outside If TRUE geosymbols will plot outside the specified limits set
#' by geoplot(). If FALSE, which is default, outside points will be skipped.
#' @param minsym Minimum symbol, default is "<", meaning that if levels = c(1,
#' 2), labels will be presented as < 1, 1-2, 2 <, but if minsym = " " labels
#' will be presented as 1, 1-2, 2. See also labels.resolution.
#' @param boundcheck If boundcheck != 0 those points which are out of bounds
#' are returned to the user, if boundcheck = 2 the points are also not plotted,
#' default is boundcheck = 0.
#' @param na.rm If true NA's are removed, default is true.
#' @param label.resolution the resolution (precision) of the label numbers,
#' default is 0, meaning that if levels = c(1, 2), labels will be presented <
#' 1, 1-2, 2 <, if label.resolution = 0.1 labels will be presented < 1, 1.1-2,
#' 2.1<, see also minsym.  If label.resolution = "none", the labels will
#' present the lowest number of the interval with each color.
#' @param characters A boolean variable determing whether characters are to be
#' plotted.
#' @param pch Type of symbols for each level.
#' @param marks Type of symbols for each level. The difference between marks
#' and pch is becae when making points on a plot the user can either give pch =
#' 17 or pch = "A". The former type is called marks here but the latter pch.
#' Marks have to be given for each level and set to -1 where pch is to be used.
#' The length of the vector
#' @param charcol The color of the charchters, default is the same as col.
#' @param open.circles Should open.circles be plotted. Default FALSE. NEEDS
#' CHECKING.
#' @param col.names Column names with positions. Default \code{lat, lon}.
#' @param border Should border be plotted? Default FALSE. NEEDS CHECKING.
#' @param bordercol Color of border. Default 0. NEEDS CHECKING.
#' @return No values are returned
#' @seealso \code{\link{geoplot}}, \code{\link{geopolygon}},
#' \code{\link{points}}, \code{\link{geotext}}, \code{\link{geopoints}},
#' \code{\link{geocontour.fill}}, \code{\link{geogrid}},
#' \code{\link{geocontour}}.
#' @keywords aplot
#' @examples
#' 
#'  \dontrun{     # lodna.2 composes of echo measurements for capelin on the
#'       # norther- and easternshores of Iceland. [lat, lon, z]
#' 
#'       # Show points.
#' 
#'       geoplot(lodna.2, type = "l", grid = FALSE)   # Begin by plotting Iceland.
#'       geopoints(lodna.2, pch = "*", col = 150)     # See where the points are.
#' 
#'       ####################################
#'       # Example 1, color parbars plot.   #
#'       ####################################
#' 
#'       geoplot(lodna.2, type = "l", grid = TRUE)     # Begin by plotting Iceland.
#'       levels = c(0, 20, 50, 100, 500, 1000)
#' 
#'       geosymbols(lodna.2, z = lodna.2$z, colplot = TRUE, colors = seven.col,
#'                  levels = levels, parbars = 0.05, colors = seven.col,
#'                  label.location = "locator")
#' 
#'       # "locator" click twice on the map where you want the contour index.
#'       # Indicate firstly the upper left corner position then lower right.
#' 
#'       #######################################
#'       # Example 2, black/white perbars plot.#
#'       #######################################
#' 
#'       geoplot(lodna.2, type = "l", grid = FALSE)
#' 
#'       geosymbols(lodna.2, z = lodna.2$z, perbars = 0.1)
#' 
#'       # Bars perpendicular to measurement direction
#' 
#'       #######################################
#'       # Example 3, Color Dots.              #
#'       #######################################
#' 
#'       # Set up data.
#'       attach("/usr/local/reikn/SplusNamskeid")
#'       i<-utbrteg$ar == 2004
#' 
#'       # Set up the plot.
#'       geoplot()
#'       levels = c(10, 100, 500)
#'       colors = c(13, 55, 111, 153)
#'       labloc<-list(lat = c(63.95, 65.4), lon = c(-19.8, -17.3))
#' 
#'       geosymbols(utbrteg[i, ], z = utbrteg[i, "torskur.kg"], circles = 0.05,
#'                  sqrt = TRUE, colplot = TRUE, levels = levels, colors = colors,
#'                  label.location = labloc)
#' 
#' 
#'       #######################################
#'       # Example 4, Rings around points.     #
#'       #######################################
#' 
#'       # Having done the set up data and plot in Example 3.
#' 
#'       geoplot(utbrteg$lat, utbrteg$lon, pch = ".")
#'       geosymbols(utbrteg[i, ], z = utbrteg[i, "torskur.kg"], circles = 0.2,
#'                  sqrt = TRUE, label.location = labloc)
#' 
#'       # Circles can be replaced with squares, rectangles, vbars, hbars or
#'       # perbars or more than one used simultanuously.
#' }
#' @export geosymbols
geosymbols <-
function(lat, lon = 0, z, levels = NULL, reflevels = NULL, labels.only = FALSE,
	cex = 0.6, chs = 0.8, z1 = 0, circles = 0, squares = 0, rectangles = c(
	0, 0), vbars = 0, hbars = 0, perbars = 0, parbars = 0, sqrt = FALSE, col = 
	1, maxn = 0, colplot = FALSE, nlevels = 10, colors = 0, n = 25, maxcol = 
	155, only.positive = FALSE, digits = 0, white = FALSE, lwd = 1, label.location
	 = NULL, labels = 1, fill.circles = FALSE, density = 0, angle = 45, rotate
	 = 0, outside = FALSE, minsym = "<", boundcheck = 0, na.rm = TRUE, 
	label.resolution = 0, characters = FALSE, pch, marks, charcol = 0, 
	open.circles = FALSE, col.names = c("lat", "lon"), border = FALSE, bordercol = 
	0)
{
	geopar <- getOption("geopar")
	options(warn = -1)
	if(!is.null(label.location))
		if(is.list(label.location))
			label.location <- as.data.frame(label.location)
	if(is.data.frame(lat)) {
		i <- match(col.names, names(lat))
		data <- data.frame(lat = lat[, i[1]], lon = lat[, i[2]])
	}
	else {
		data <- data.frame(lat = lat, lon = lon)
	}
	if(na.rm) {
		# delete na.
		ind <- c(1:length(data$lat))
		ind <- ind[is.na(data$lat) | is.na(data$lon)]
		if(length(ind) > 0) {
			data$lat <- data$lat[ - ind]
			data$lon <- data$lon[ - ind]
			z <- z[ - ind]
		}
	}
	ind <- c(1:length(data$lat))
	ind <- ind[is.na(z)]
	if(length(ind) > 0) {
		data$lat[ind] <- NA
		data$lon[ind] <- NA
		z[ind] <- mean(z, na.rm = TRUE)
	}
	if(fill.circles)
		colplot <- TRUE
	if(open.circles)
		colplot <- TRUE
	if(density > 0)
		colplot <- TRUE
	if(maxn == 0)
		maxn <- max(abs(z))
	if(only.positive) {
		ind <- c(1:length(z))
		ind <- ind[z < 0]
		z[ind] <- 0
	}
	if(boundcheck != 0) {
		dataprj <- Proj(data$lat, data$lon)
		ind <- c(1:length(data$lat))
		ind <- ind[dataprj$x < geopar$limx[1] | dataprj$x > geopar$
			limx[2] | dataprj$y < geopar$limy[1] | dataprj$y > 
			geopar$limy[2]]
		if(length(ind) > 0) {
			ind1 <- paste(ind, collapse = ",")
			print(paste("points", ind1, "out of bounds"))
		}
		if(boundcheck == 2) {
			if(length(ind) > 0) {
				data$lat <- data$lat[ - ind]
				data$lon <- data$lon[ - ind]
				z <- z[ - ind]
			}
		}
	}
	oldpar <- selectedpar()
	par(geopar$gpar)
	on.exit(par(oldpar))
	if(outside)
		par(xpd = TRUE)
	else par(xpd = FALSE)
	if(colplot) {
		if(labels.only) {
			if(cex != 0)
				par(cex = cex)
			colsymbol(data$lat, data$lon, z, circles, squares,
				rectangles, hbars, vbars, perbars, parbars,
				levels, nlevels, colors, white, n, maxcol,
				digits, label.location, labels, fill.circles,
				density, angle, rotate, minsym, 
				label.resolution, col, labels.only = TRUE, 
				open.circles = open.circles, lwd = lwd, border
				 = border, bordercol = bordercol)
		}
		else {
			if(cex != 0)
				par(cex = cex)
			colsymbol(data$lat, data$lon, z, circles, squares,
				rectangles, hbars, vbars, perbars, parbars,
				levels, nlevels, colors, white, n, maxcol,
				digits, label.location, labels, fill.circles,
				density, angle, rotate, minsym, 
				label.resolution, col, open.circles = 
				open.circles, lwd = lwd, border = border, 
				bordercol = bordercol)
		}
	}
	else {
		x <- Proj(data$lat, data$lon, geopar$scale, geopar$b0, geopar$
			b1, geopar$l1, geopar$projection)
		y <- x$y
		x <- x$x
		ein.pr.in <- (geopar$limy[2] - geopar$limy[1])/geopar$gpar$
			pin[2]
		if(!is.null(label.location)) {
			if(label.location == "locator" || label.location == 0)
				label.location <- geolocator(n = 2)
			limits <- Proj(label.location)
			xlim <- limits$x
			ylim <- limits$y
			if(xlim[1] > xlim[2]) {
				temp <- xlim[1]
				xlim[1] <- xlim[2]
				xlim[2] <- temp
			}
			if(ylim[1] > ylim[2]) {
				temp <- ylim[1]
				ylim[1] <- ylim[2]
				ylim[2] <- temp
			}
			if(is.null(levels)) {
				rg <- range(z)
				lrg <- rg[2] - rg[1]
				levels <- signif(seq(rg[1] + lrg/10, rg[2] -
					lrg/10, length = nlevels), 2)
			}
			lbox <- length(levels)
			boxy <- c(1:lbox)
			boxy <-  - boxy/lbox + 1
			boxy1 <- boxy + 1/(1.2 * lbox)
			yloc <- (boxy + boxy1)/2
			xloc <- matrix(0.85, length(yloc))
			par(adj = 0)
			textx <- as.character(levels)
			boxx <- c(matrix(0.1, 1, length(boxy)))
			boxx <- xlim[1] + abs((xlim[2] - xlim[1])) * boxx
			xloc <- xlim[1] + abs((xlim[2] - xlim[1])) * xloc
			yloc <- ylim[1] + abs((ylim[2] - ylim[1])) * yloc
			boxy <- ylim[1] + (ylim[2] - ylim[1]) * boxy
			ll <- (ylim[2] - ylim[1]) * 0.05
			if(circles != 0 | squares != 0 | hbars != 0 | vbars !=
				0 | perbars != 0)
				text(boxx, boxy + ll, textx, cex = chs)
		}
		if(circles != 0) {
			# plot circles.
			rg <- range(circles)
			rglen <- rg[2] - rg[1]
			lev <- seq(rg[1] + rglen/10, rg[2] - rglen/10, length
				 = 5)
			if((circles > 100) | (circles < 0))
				circles <- 0.2
			#default value.  
			circles <- ein.pr.in * circles
			# size in units
			if(sqrt) {
				if(!labels.only)
					symbols(x, y, circles = circles * sqrt(
						abs(z)/maxn), inches = FALSE, add
						 = TRUE, fg = col, lwd = lwd)
				if(!is.null(label.location))
					symbols(c(xloc), c(yloc), circles = 
						circles * sqrt(abs(levels)/
						maxn), add = TRUE, inches = FALSE,
						lwd = lwd, fg = col)
			}
			else {
				if(!labels.only)
					symbols(x, y, circles = circles * (
						abs(z)/maxn), add = TRUE, inches
						 = FALSE, fg = col, lwd = lwd)
				if(!is.null(label.location))
					symbols(c(xloc), c(yloc), circles = (
						circles * abs(levels))/maxn,
						add = TRUE, inches = FALSE, lwd = lwd, fg = col)
			}
		}
		if(squares != 0) {
			#plot squares.
			if((squares > 100) | (squares < 0)) squares <- 0.2
			#default value.  
			squares <- ein.pr.in * squares
			# size in units  
			if(sqrt) {
				if(!labels.only)
					symbols(x, y, squares = squares * sqrt(
						abs(z)/maxn), add = TRUE, inches
						 = FALSE, fg = col, lwd = lwd)
				symbols(c(xloc), c(yloc), squares = squares *
					sqrt(abs(levels)/maxn), add = TRUE, inches
					 = FALSE, lwd = lwd, fg = col)
			}
			else {
				if(!labels.only)
					symbols(x, y, squares = squares * (
						abs(z)/maxn), add = TRUE, inches
						 = FALSE, fg = col, lwd = lwd)
				symbols(c(xloc), c(yloc), squares = (squares *
					abs(levels))/maxn, add = TRUE, inches = FALSE, fg = col,
					lwd = lwd)
			}
		}
		if((rectangles[1] != 0) | (rectangles[2] != 0)) {
			# plot rectangles
			if((rectangles[1] > 100) | (rectangles[1] < 0)) 
					rectangles[1] <- 0.2
			# thickness
			if((rectangles[2] > 100) | (rectangles[2] < 0)) 
					rectangles[2] <- 0.2
			# length
			rectangles[1] <- ein.pr.in * rectangles[1]
			# size in units  
			rectangles[2] <- ein.pr.in * rectangles[2]
			# size in units  
			m <- matrix(rectangles[1], length(z), 2)
			if(sqrt)
				m[, 1] <- rectangles[1] * sqrt(abs(z)/maxn)
			else m[, 1] <- (rectangles[1] * abs(z))/maxn
			if(length(z1) > 1) {
				if(sqrt)
					m[, 2] <- rectangles[2] * sqrt(abs(
						z1)/max(abs(z1)))
				else m[, 2] <- (rectangles[2] * abs(z1))/max(
						abs(z1))
			}
			symbols(x, y, rectangles = m, add = TRUE, inches = FALSE,
				fg = col, lwd = lwd)
		}
		if(vbars != 0) {
			# plot vertical bars
			if(vbars > 100) vbars <- 0.4
			mx <- matrix(NA, 3, length(x))
			my <- mx
			mx[1,  ] <- x
			my[1,  ] <- y
			mx[2,  ] <- x
			#      mlocx<- matrix(NA, 3, length(levels)); mlocy<-mlocx
			#      mlocx[1, ]<-c(xloc) ; mlocy[1, ]<-c(yloc); mlocx[2, ]<-c(xloc)
			r <- ein.pr.in * vbars
			# size in units  
			if(sqrt) {
				my[2,  ] <- my[1,  ] + r * sqrt(abs(z)/maxn)
			}
			else {
				my[2,  ] <- my[1,  ] + (r * abs(z))/maxn
			}
			if(!labels.only)
				lines(mx, my, col = col, lwd = lwd)
		}
		if(hbars != 0) {
			# plot horizontal bars
			if(hbars > 100) hbars <- 0.4
			mx <- matrix(NA, 3, length(x))
			my <- mx
			mx[1,  ] <- x
			my[1,  ] <- y
			my[2,  ] <- y
			#      mlocx<- matrix(NA, 3, length(levels)); mlocy<-mlocx
			#      mlocx[1, ]<-c(xloc) ; mlocy[1, ]<-c(yloc); mlocy[2, ]<-c(yloc)
			r <- ein.pr.in * hbars
			# size in units
			if(sqrt) {
				mx[2,  ] <- mx[1,  ] + r * sqrt(abs(z)/maxn)
			}
			else {
				mx[2,  ] <- mx[1,  ] + (r * abs(z))/maxn
			}
			if(!labels.only)
				lines(mx, my, col = col, lwd = lwd)
		}
		if(perbars != 0) {
			# plot bars perpendicular to cruiselines
			if(perbars > 100) perbars <- 0.4
			mx <- matrix(NA, 3, length(x))
			my <- mx
			mx[1,  ] <- x
			my[1,  ] <- y
			#      mlocx<- matrix(NA, 3, length(levels)); mlocy<-mlocx
			#      mlocx[1, ]<-c(xloc) ; mlocy[1, ]<-c(yloc); mlocy[2, ]<-c(yloc)
			r <- ein.pr.in * perbars
			# size in units  
			dx <- c(1:length(x))
			dx[1] <- x[2] - x[1]
			dx[2:(length(x) - 1)] <- x[3:(length(x))] - x[1:(length(
				x) - 2)]
			dx[length(x)] <- x[length(x)] - x[length(x) - 1]
			dy <- c(1:length(y))
			dy[1] <- y[2] - y[1]
			dy[2:(length(y) - 1)] <- y[3:length(y)] - y[1:(length(
				y) - 2)]
			dy[length(y)] <- y[length(x)] - y[length(y) - 1]
			dxy <- sqrt(dx * dx + dy * dy)
			dx <- dx/dxy
			dy <- dy/dxy
			if(sqrt)
				mx[2,  ] <- mx[1,  ] - dy * r * sqrt(abs(z)/
					maxn)
			else mx[2,  ] <- mx[1,  ] - (dy * r * abs(z))/maxn
			if(sqrt)
				my[2,  ] <- my[1,  ] + dx * r * sqrt(abs(z)/
					maxn)
			else my[2,  ] <- my[1,  ] + (dx * r * abs(z))/maxn
			#      if(sqrt)  mlocx[2, ]<-mlocx[1, ]+r*sqrt(abs(levels)/maxn)
			#      else  mlocx[2, ]<-mlocx[1, ]+r*abs(z)/maxn
			if(!labels.only) lines(mx, my, col = col)
		}
	}
	par(oldpar)
	if(characters) {
		if(missing(marks))
			marks <- rep(-1, length(pch))
		if(missing(pch))
			pch <- rep(" ", length(marks))
		if(!is.numeric(levels)) {
			if(!is.numeric(z))
				ind <- match(z, levels)
			else {
				if(is.null(reflevels)) {
					print("Error")
					return()
				}
				ind <- match(z, reflevels)
			}
			n <- length(levels)
		}
		else {
			if(charcol != 0)
				col <- charcol
			n <- length(levels) + 1
			levels <- c(-1000000., levels, 1000000.)
			ind <- cut(z, levels, labels = FALSE)
		}
		if(length(col) == 1)
			col <- rep(col, n)
		if(length(cex) == 1)
			cex <- rep(cex, n)
		if(!labels.only) {
			for(i in 1:n) {
				tmp <- data[ind == i,  ]
				if(nrow(tmp) > 0) {
					if(marks[i] < 0)
						geopoints(tmp, pch = pch[i],
							cex = cex[i], col = col[
							i])
					else geopoints(tmp, pch = marks[i],
							cex = cex[i], col = col[
							i])
				}
			}
		}
		if(!is.null(label.location)) {
			########
			if(!is.list(label.location)) if(label.location == 
					"locator")
					label.location <- geolocator(n = 2)
			oldpar <- selectedpar()
			on.exit(par(oldpar))
			par(geopar$gpar)
			paint.window(label.location)
			label.location <- Proj(label.location)
			## if(is.numeric(levels))
			## 	Pointlabel(levels[2:(length(levels) - 1)],
			## 		digits, label.location$x, 
			## 		label.location$y, minsym, 
			## 		label.resolution, marks, pch, col,
			## 		cex, chs)
			## else Charlabel(levels, label.location$x, label.location$
			## 		y, label, marks, pch, col, cex, chs)
		}
	}
	options(warn = 0)
	return(invisible())
}

