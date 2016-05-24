#' Function that fills rectangles or sub-rectangles according to the amount of
#' fish calculated in them.
#' 
#' Function that fills rectangles or sub-rectangles according to the amount of
#' fish calculated in them, or some other quantity.
#' 
#' Geoplot has to be called before to set the parameter geopar.  If geoplot is
#' called with the parameter cont = TRUE labels are plotted on the left hand
#' side of the plot. If a color plot is to be made litaps has to be called
#' before.
#' 
#' @param reitur Number of square.  If number of square >1000 the last digit
#' means subsquare.  For example 5604 means square 560, subsquare 4.
#' @param z Values for each square and subsquare.
#' @param smareitur Number of subsquare if they are used.  Default is that
#' subsquare is not used.
#' @param levels Values between colors or different density of lines.  Default
#' is to let the program find it. A small example shows how levels are given.
#' <s-example>
#' 
#' levels <- c(0, 100, 1000, 2000, 4000, 5000, 8000)
#' 
#' </s-example> Then the groups are <100, 100-1000, 2000-4000, 4000-5000 &
#' >5000, i.e.  the program more or less ignores the first and the last values
#' given.
#' @param colors Number of color or shading pattern for each level.  Colours
#' can go from 0 thru 155.  They go from white = 0, black = 1, blue = 2, green,
#' yello, red = 155.  Density goes from 10 thru 130.  Everything above 130 is
#' completely black.  If colors is not given the program determines it from the
#' data.  Default maximum density is set to 70 but everything above that is
#' very black.
#' @param density If density is 0 coloured plot is made else different desity
#' of lines is used.
#' @param maxcol Maximum color ( or density) Default value is 155 for colors 70
#' for shading.  All density of lines > 130 are black.
#' @param nlevels Number of levels.  Default value is 6.
#' @param white If white = TRUE the lowest color or density is zero.  Default
#' value is white.
#' @param border If border of polygon is to be drawn.  Default value is FALSE.
#' @param angle Angle of lines used for hatching.  Default value is 45 degrees.
#' @param rotate Rotation of hatcing lines between each level.  Default value
#' is zero or no rotation but rotate = 90 is very useful.
#' @param digits Number of digits in labels.  Default value is zero.
#' @param cex Charcter expansion of digits in labels.  Default value is the
#' current size when the program is called.
#' @param label.location location of labels, default is no labels are plotted.
#' @param labels.only if true only labels are plotted.
#' @param col the color of lines, default is 1 (black).
#' @param outside if true squares outside plot area are plotted, default is
#' true.
#' @param mincol the color number of the first color to be used, colors are
#' taken from the range mincol to maxcol.
#' @return No values returned.
#' @section Side Effects: No side effects.
#' @seealso \code{\link{geoplot}}.
#' @examples
#' 
#'        geoplot(reitur = TRUE)
#'        reitur <- c(560, 560, 560, 560, 561, 561, 561, 561)
#'        smareitur <-c(1, 2, 3, 4, 1, 2, 3, 4)
#'        lab.loc <- list(lat = c(63.95, 65.4), lon = c(-19.8, -17.3))
#'        z <- c(100, 200, 300, 222, 345, 453, 232, 212)
#'        geopolygon(island); geolines(island)
#'        reitaplott(reitur, smareitur, z, rotate = 90, label.location = lab.loc)
#'        # All subsuares in square 560 and 561.
#' 
#' @export reitaplott
reitaplott <-
function(reitur, smareitur = NULL, z, levels = NULL, colors = NULL, density = 1,
	maxcol = 155, nlevels = 6, white = FALSE, border = FALSE, angle = 45, rotate = 
	0, digits = 0, cex = 0, label.location = 0, labels.only = FALSE, col = 1,
	outside = TRUE, mincol = NULL)
{
	geopar <- getOption("geopar")
	oldpar <- selectedpar()
	par(geopar$gpar)
	on.exit(par(oldpar))
	brd <- list(x = c(geopar$limx[1], geopar$limx[2], geopar$limx[2], 
		geopar$limx[1], geopar$limx[1]), y = c(geopar$limy[1], geopar$
		limy[1], geopar$limy[2], geopar$limy[2], geopar$limy[1]))
	if(geopar$cont)
		par(plt = geopar$contlines)
	z <- z + 1e-08
	brd$x <- brd$x + 1e-06
	brd$y <- brd$y + 1e-07
	if(abs(reitur[1]) > 999 && is.null(smareitur)) {
		smareitur <- abs(reitur) - floor(abs(reitur)/10) * 10
		reitur <- sign(reitur) * floor(abs(reitur)/10)
	}
	or.angle <- angle
	dlat <- c(0.125, 0.125, -0.125, -0.125)
	dlon <- c(-0.25, 0.25, -0.25, 0.25)
	lat <- r2d(reitur)
	lon <- lat$lon
	lat <- lat$lat
	if(!is.null(smareitur)) {
		lat <- lat + dlat[smareitur]
		lon <- lon + dlon[smareitur]
		dlon <- c(-0.25, 0.25, 0.25, -0.25, -0.25)
		dlat <- c(-0.125, -0.125, 0.125, 0.125, -0.125)
	}
	else {
		dlon <- c(-0.5, 0.5, 0.5, -0.5, -0.5)
		dlat <- c(-0.25, -0.25, 0.25, 0.25, -0.25)
	}
	nlat <- length(lat)
	nlon <- length(lon)
	lat <- t(matrix(lat, nlat, 5))
	lon <- t(matrix(lon, nlon, 5))
	dlon <- matrix(dlon, 5, nlat)
	dlat <- matrix(dlat, 5, nlat)
	lat <- lat + dlat
	lon <- lon + dlon
	n.a <- rep(NA, ncol(lat))
	lat <- rbind(lat, n.a)
	lon <- rbind(lon, n.a)
	# 	Set levels and color.  	
	z <- z + 1e-07
	# because of zeroes.  
	if(is.null(levels)) {
		if(nlevels == 0)
			nlevels <- 10
		levels <- pretty(range(z, na.rm = TRUE), nlevels)
		levels <- levels[2:(length(levels) - 1)]
	}
	ncont <- length(levels)
	if(is.null(colors)) {
		if(density > 0 && is.null(mincol))
			mincol <- 8
		if(density == 0 && is.null(mincol))
			mincol <- 2
		if(density > 0 && is.null(maxcol))
			maxcol <- 70
		if(density == 0 && is.null(maxcol))
			maxcol <- 155
		if(white) {
			# lowest values white.  
			colors <- c(1:(ncont))
			colors <- floor(mincol + ((colors - 1) * (maxcol - 
				mincol))/(length(colors) - 1))
			colors <- c(0, colors)
		}
		else {
			colors <- c(1:(ncont + 1))
			colors <- floor(mincol + ((colors - 1) * (maxcol - 
				mincol))/(length(colors) - 1))
		}
	}
	print(paste("calculated colors", paste(colors, collapse = ",")))
	print(paste("calculated levels", paste(levels, collapse = ",")))
	levels.1 <- levels
	colors.1 <- colors
	m <- max(z[!is.na(z)])
	if(!is.null(levels)) {
		i <- c(1:length(levels))
		i <- i[levels > max(z[!is.na(z)])]
		if(length(i) > 0) {
			levels <- levels[ - i]
			if(length(colors) > 1) {
				i <- i + 1
				colors <- colors[ - i]
			}
		}
	}
	levels <- c(levels, max(c(max(abs(levels)) * 1.1, max(z) * 1.1)))
	# change.
	levels <- c(min(c(min(z[z != -99999]) - 1, levels[1] - 1)), levels)
	if(levels[2] - levels[1] < 1)
		levels[1] <- levels[2] - 1
	print("adjusted contour are")
	print(round(levels, 3))
	print("adjusted colors are")
	print(round(colors, 3))
	#	Define color for each point.  
	ind <- cut(z, levels)
	ind1 <- ind
	ind <- colors[ind]
	# number of color. 
	x <- Proj(lat, lon, geopar$scale, geopar$b0, geopar$b1, geopar$l1,
		geopar$projection)
	if(density > 0 || !outside) {
		rot <- rep(0, length(colors))
		rot[1] <- angle
		for(i in 2:length(colors))
			rot[i] <- rot[i - 1] + rotate
		rot <- rot[ind1]
		if(!labels.only) {
			for(i in 1:length(ind)) {
				if(geopar$projection == "Lambert")
					x1 <- fill.points(x$x[1:5, i], x$y[
						1:5, i], nx = 10, option = 2)
				else x1 <- list(x = x$x[1:5, i], y = x$y[1:
						5, i])
				if(!outside)
					x1 <- findcut(x1, brd)
				if(length(x1$x) > 3) {
					if(density > 0) {
						if(ind[i] > 100)
							polygon(x1$x, x1$y,
								col = 1, border
								 = border)
						else polygon(x1$x, x1$y, 
								density = ind[
								i], border = 
								border, angle
								 = rot[i])
					}
					else polygon(x1$x, x1$y, col = ind[
							i], border = border)
					if(border) {
						if(!outside)
							x1 <- findline(x1,
								brd)
						if(length(x1$x) != 0)
							lines(x1$x, x1$y, col
								 = col)
					}
				}
			}
		}
		#       Add  labels around plot
		if(geopar$cont) {
			par(plt = geopar$contlab)
			par(new = TRUE)
			plot(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0), type = "l",
				axes = FALSE, xlab = " ", ylab = " ")
			shading1(levels.1, digits, colors, angle = or.angle,
				rotate = rotate, cex = cex, fill = geopar$
				cont)
		}
	}
	else {
		# use colors 
		if(!labels.only) {
			if(geopar$projection == "Lambert")
				x <- fill.points(x$x, x$y, 10, option = 2)
			polygon(x$x, x$y, col = ind, border = border)
			if(border == TRUE)
				lines(x$x, x$y)
		}
		#       Add  labels around plot
		if(geopar$cont) {
			par(plt = geopar$contlab)
			par(new = TRUE)
			plot(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0), type = "l",
				axes = FALSE, xlab = " ", ylab = " ")
			if(cex != 0)
				par(cex = cex)
			labels1(levels.1, digits, colors.1, fill = geopar$
				cont)
		}
	}
	# 	Add  labels around plot 
	if(length(label.location) == 1) if(label.location == "locator")
			label.location <- geolocator(n = 2)
	if(length(label.location) > 1) {
		#label located somewhere in drawing
		label.location <- Proj(label.location$lat, label.location$
			lon, geopar$scale, geopar$b0, geopar$b1, geopar$l1,
			geopar$projection)
		paint.window(label.location, border = TRUE)
		if(density > 0)
			shading1(levels.1, digits, colors.1, angle = or.angle,
				rotate = rotate, cex = cex, xlim = 
				label.location$x, ylim = label.location$y)
		else {
			labels1(levels.1, digits, colors.1, xlim = 
				label.location$x, ylim = label.location$y)
		}
	}
	par(oldpar)
	return(invisible())
}

