#' geocontour.fill plots colored or black and white contours on a graph made by
#' geoplot.
#' 
#' The program accepts data on a rectangular grid.  Irregular data can be
#' interpolated on the grid in a number of ways, using kriging, the Splus
#' function interp, etc.  NA's are allowed in the matrix but they are changed
#' to 0 or the average value of z early in the program. The geo package
#' contains two different contour plot programs, geocontour and
#' geocontour.fill.  Geocontour draws contour lines and geocontour.fill fills
#' in between the lines.  Hardcopy of the plot can be done on a color
#' postscript printer or on a black and white postscript printer.  In all cases
#' the "postscript" driver included with Splus has to be used.  There is
#' another postscript driver included with Splus called "pscript".  That driver
#' is used when a plot is on the screen and the print button is pressed.  That
#' can not be done when the postscript driver is used as described here after.
#' 
#' 
#' @param grd List with components \$lat and \$lon (or \$x, \$y) defining the
#' grid values.  Can be made for example by lat <- list(lat = seq(60, 65,
#' length = 100), lon = seq(-30, - 10, length = 200)).  Also lat can be the
#' outcome of the program grid and then the components lat\$grpt\$lat,
#' lat\$grpt\$lon and lat\$reg are used.
#' @param z Matrix or vector of values.  Length(z) = nlat*nlon except when lat
#' is a list with components lat\$grpt and lat\$xgr.  Then the length of z is
#' the same as the length of lat\$xgr\$lat and lat\$xgr\$lon.  This is the case
#' if the output of the grid program is used as input lat.
#' @param levels Values at contourlines.  Default value is zero.  If levels is
#' zero the program determines the contourlines from the data.  If levels is of
#' length 3 with levels = c(0, 1, 2) then there are 4 groups each with special
#' color i.e. <0, 0-1, 1-2 and > 2.
#' @param nlevels Number of contourlines.  Used if the program has to determine
#' the contourlines.  Default value is 10.  nlevels = 10 does not always mean
#' 10 due to characteristic of the pretty command.
#' @param cex Character size expansion.  Size of letters in labels.  Default
#' value is 0.7.
#' @param digits Number of digits in the labels.  Default value is 1.
#' @param col Color number for the contourlines.  A vector one longer than the
#' vector levels Default is blue-green- yellow-red from lowest to the highest
#' values.  The program assumes certain setup of the splus colors that is
#' described later. Colors should not be specified explicitly except the
#' contour lines are.
#' @param working.space Size of working space.  The program determines it from
#' data and prints on the screen.  If mysterious errors occur then it is likely
#' that the program has not reserved enough working space.  The program writes
#' the used work space on the screen for use in similar situations. (saves
#' time.)
#' @param labels Type of labels, either 1 or 2 <s-example>
#' 
#' labels = 1 means type of label used with few contourlines <20.  labels = 2
#' means type of label used with many contourlines >20.  can also use labels =
#' 3 means only labels, no contour plot.
#' 
#' </s-example> Labels can be inserted in two ways, by calling geoplot with
#' cont = TRUE or by specifying label.location.  In the first case the left
#' part of the plot is reserved for labels while in the latter case the label
#' is put in a place specified by the user.  The latter method is recommended
#' in most cases.
#' @param ratio Factor used to avoid numerical problems.  Default value 100.
#' Lower values decrease numerical problems but can introduce bias.
#' @param only.positive Logical value.  If TRUE then negative values are not
#' allowed else negative values are set to zero.  Default value is FALSE.
#' @param fill Determines whether NA's should be replaced with zeros (fill = 1)
#' or mean(z) (fill = 2).  Default value is fill = 1.  Used when lat is a list
#' with components lat\$xgr\$lat, lat\$xgr\$lon, lat\$grpt\$lat and
#' lat\$grpt\$lon.  Only used in special cases.
#' @param maxcol Number of colors used (excluding #0).  Default value 155.
#' @param white If true the the first class is represented with white (color
#' 0). Default value is FALSE but white = TRUE is also often used.
#' @param label.location List with components \$lat and \$lon.  (or \$x, \$y)
#' Gives the lower left and upper right corner of the box where the labels are
#' put.  Default value is 0 that means no labels are put on the drawing (except
#' when geopar\$cont = TRUE).  l1 is best given by geolocator or directly by
#' specifying label.location = "locator".
#' @param labels.only If labels.only is true no contours are drawn but only
#' labels.  Default value is false.  Order of the commands could be:
#' <s-example> > geoplot(deg, plot = FALSE) > geocontour.fill(lat, z = z,
#' levels = lev, reg = area) # Draw contours.  > geocontour(lat, z = z, levels
#' = lev, colors = FALSE, reg = area) # Only used with black and white printers
#' to make distinction # between different levels clearer.  (not used with
#' color # printers).  > geoplot(deg, new = TRUE) # Refresh gridlines.  >
#' geocontour.fill(lat, z = z, levels = lev, labels.only = TRUE, label.location
#' = l1) # Add labels.  </s-example>
#' @param bordercheck if true the program checks if plot is outside border and
#' does not plot what is outside border.
#' @param maxn a parameter determing the resolution of the plot, unimportant.
#' @param bcrat bordercheck ratio, how much outside the border we will allow to
#' be plotted, default is bcrat = 0.05 meaning that we will allow the plot to
#' go 5\% off the border.
#' @param limits To be described.
#' @param col.names the names of the columns in grid, the first argument will
#' be plotted on the x-axes and the second on the y-axes.  Default is col.names
#' = c("lat", "lon").
#' @param minsym minimum symbol, default is "<", meaning that if levels = c(1,
#' 2), labels will be presented as < 1, 1-2, 2 <, but if minsym = " " labels
#' will be presented as 1, 1-2, 2.  See also labels.resolution.
#' @param label.resolution the resolution (precision) of the label numbers,
#' default is 0, meaning that if levels = c(1, 2), labels will be presented <
#' 1, 1-2, 2 <, if label.resolution = 0.1 labels will be presented < 1, 1.1-2,
#' 2.1<, see also minsym.  If label.resolution = "none", the labels will
#' present the lowest number of the interval with each color.
#' @param labtxt To be described.
#' @param boxcol Colour of box around labels (legend?).
#' @param first.color.trans To be described.
#' @param mai \code{par} argument 'margin in inches'?
#' @param leftrat To be described.
#' @param labbox should a box be drawn around labels (legend?).
#' @param csi Size of character.  This parameter can not be set in R but for
#' compatibility with old Splus scripts the parameter cex is readjusted by cex
#' = cex*csi/0.12.  Use of this parameter is not recommended.  Default value is
#' NULL i.e not used.
#' @section Details: <s-example>
#' 
#' The program is based on making triangles out of a matrix of data nx * ny.
#' The number of triangles is (nx-1)*(ny- 1)*4 and linear interpolation is used
#' inside each triangle.  The factor 4 comes from the fact that the number of
#' points is doubled by interpolation. </s-example> <s-example>
#' 
#' Most of the calculations in the program are done by a program written in C.
#' The Splus part of the program prepares the data for the C program, reserves
#' memory and the user interface is in Splus.  The C program is loaded
#' automatically. </s-example> <s-example>
#' 
#' The program is currently based on the following color setup in Splus.
#' </s-example> <s-example>
#' 
#' splus*color : white black blue 50 green 50 yellow 50 red </s-example>
#' <s-example>
#' 
#' This gives 155 colors with color#0 white, #1 black, #2 blue, #53 green, #104
#' yellow, and #155 red.  The parameter maxcol is set to 155 based on this
#' setup but it can be changed if another setup is used.  The vector postcol
#' stores the mapping from the colors on the terminal to color postscript based
#' on this setup. </s-example> <s-example>
#' 
#' If the openlook() or motif() windowmanager is used in Splus the colors can
#' be changed inside Splus.  In openlook it is best to make a colorscheme
#' called geoplot: black blue 50 green 50 yellow 50 red.  The background color,
#' i.e. white is not in the definition here. </s-example> <s-example>
#' 
#' The maximum number of contourlines that can be used is currently 60.
#' </s-example> <s-example>
#' 
#' The program requires that geoplot is called before to set up the drawing.
#' geoplot is called by the parameter cont = TRUE. This is to reserve space for
#' labels on the drawing. </s-example> <s-example>
#' 
#' If the drawing is written to a file geoplot should be called with plot =
#' FALSE.  Then it only sets up the drawing but does not plot anything so the
#' size of the file will be reduced. </s-example> <s-example>
#' 
#' If the labels do not fit the relative space taken by labels and picture can
#' be changed by the parameter lcont in geoplot.  After geocontour.fill is
#' called geoplot is called again with new = TRUE to get all kind of lines back
#' but they were painted over by the program. </s-example> <s-example>
#' 
#' The program treats borders in a special way.  It begins by making
#' contour-lines over the hole area as if the borders did not exist.  When that
#' is finished the area outside the borders is painted white. </s-example>
#' <s-example>
#' 
#' If good picture is needed geocontour should be used between levels to get
#' sharper pictures. </s-example> <s-example>
#' 
#' If the matrix z is not full the program should be called by a list
#' lat\$xgr\$lat, lat\$xgr\$lon, lat\$grpt\$lat and lat\$grpt\$lon.  Then the
#' length of the vectors z, lat\$xgr\$lat and lat\$xgr\$lon is the same.
#' lat\$grpt\$lat and lat\$grpt\$lon is on the other hand the coordinates of
#' the rows and columns of the matrix. The program grid makes a list with
#' components with these names. </s-example> <s-example>
#' 
#' The functions geolines, geopoint, geopolygon, geotext & geosymbols can be
#' used to add things to the contourplot. </s-example>
#' @seealso \code{\link{geoplot}}, \code{\link{geolines}},
#' \code{\link{geopolygon}}, \code{\link{geotext}}, \code{\link{geosymbols}},
#' \code{\link{geogrid}}, \code{\link{geopar}}, \code{\link{geolocator}},
#' \code{\link{geocontour}}, \code{\link{reitaplott}}, \code{\link{geodefine}}.
#' @examples
#' 
#' \dontrun{  ######################################################
#'   # Example 1.                                         #
#'   ######################################################
#' 
#'   # Need the data.frame botnv.2004 to be able to compile this
#'   # example, if not attached use: 
#'   # >attach("/usr/local/reikn/Splus5/Aflaskyrslur/Data")
#' 
#'   codgrd <-list(lat = seq(62, 68, by = 0.1), lon = seq(-30, -9, by = 0.25))
#'   # A grid is made.
#'   lab.loc<-list(lat = c(63.9, 65.6), lon = c(-20.75, -16.5)) 
#'   # Location of the label.
#' 
#'   tmp <- combine.rt(botnv.2004$lat, botnv.2004$lon,
#'                     botnv.2004$torskur, codgrd, fun = "sum", fill = TRUE)
#'   # The data is read into the grid with combine.rt.
#' 
#'   tmp$z <- tmp$z/(cos(tmp$lat*pi/180)*0.1*60*0.25*60)
#'   # Data changed.(from being in
#'   # kilos per box to kilos per square mile).
#' 
#'   vg <- list(nugget = 0.1, sill = 1, range = 50)
#'   # Parameters for the variogram
#' 
#'   z <- pointkriging(tmp$lat, tmp$lon, tmp$z, codgrd, vg,  
#'                     maxnumber = 80, maxdist = 30, set = -1) 
#'   # Dataset smoothened with pointkriging.
#' 
#'   geoplot(lat = c(63, 67.5), lon = c(-27, -11), grid = FALSE, axlabels = TRUE, type = "n")
#'   # Plot initialized
#'   level = c(160, 200, 320, 500, 700, 1000, 2000, 3000, 4000, 5000, 6000)
#'   # Levels for geocontour.fill
#' 
#'   geocontour.fill(codgrd, z, levels = level, white = TRUE      # Plot the data.
#'                 , working.space = 300000)
#' 
#'   geopolygon(island)
#'   # Contourlines inside Iceland overwritten.
#'   geolines(island)
#'   # Iceland redrawn with geolines.
#' 
#'   geocontour.fill(codgrd, z, levels = level, white = TRUE, 
#'                   label.location = lab.loc, labels.only = TRUE)
#'   # Call geocontour.fill again only to plot the labels.
#' 
#'   #########################################################
#'   # Example 2.                                            #
#'   #########################################################
#' 
#'   
#'   # Preperation for pointkriging    
#'   # th4.2002 is the data used here, dataframe [lon, lat, mat].
#'   # >attach(?????) 
#'   
#' 
#'   grd.smb <-list(lat = seq(62.8, 67.5, length = 80),   # Set up the grid.
#'                  lon = seq(-28, -10, length = 130))
#'   m.lev<-c(0, 0.1, 0.2, 0.3, 0.5)
#'   # Levels for the geocontour.fill.
#'   m.col<-c(0, 14, 59, 104, 119, 149)
#'   # Colors for the levels.
#'   lab.in.island<-list(lat = c(63.9, 65.6), lon = c(-20.75, -16.5)) 
#'   # Location of the Label
#' 
#'   vg <- list(nugget = 0.3, sill = 1, range = 50) 
#'   # Initialize variogram parameters.
#'  
#'   zfj<-pointkriging(th4.2002$lat, th4.2002$lon, z = th4.2002$mat, 
#'                   grd.smb, vg, maxnumber = 80, maxdist = 30, set = -1)
#'   # Smooth the data with pointkriging.
#' 
#'   #
#'   # Plotting                       
#'   #
#' 
#'   par(mfrow = c(1, 1),  mai = rep(0, 4))
#'   # Set up graphic parameters.
#'   geoplot(lat = c(63, 67.5), lon = c(-27, -11), grid = FALSE, axlabels = FALSE, type = "n")
#'   # Draw a background with Iceland with geoplot.
#' 
#'   geocontour.fill(zfj, levels = m.lev, col = m.col, working.space = 300000) 
#'   # Plot the data with geocontour.fill.
#' 
#'   geopolygon(gbdypif.500, col = 0, exterior = TRUE, r = 0)
#'   # Remove contours outside gbdypif.500.
#'   geoplot(lat = c(63, 67.5), lon = c(-27, -11), grid = FALSE,          # Replot.
#'           axlabels = FALSE, type = "n", new = TRUE)
#' 
#'   geopolygon(island, col = 43)
#'   # Remove contours inside Iceland and color Iceland.
#'   geocontour.fill(zfj, levels = m.lev, label.location = lab.in.island, 
#'                   labels.only = TRUE, csi = 0.1, col = m.col, working.space = 300000)
#'   # Call geocontour.fill to plot labels.
#'   geolines(island)
#'   # Redraw the lines of Iceland.
#' }
#' @export geocontour.fill
geocontour.fill <-
function(grd, z, levels = NULL, nlevels = 0, cex = 0.7, digits = 1, col = NULL,
	working.space = 0, labels = 1, ratio = 1000, only.positive = FALSE, fill = 
	0, maxcol = 155, white = FALSE, label.location = 0, labels.only = FALSE, 
	bordercheck = FALSE, maxn = 10000, bcrat = 0.05, limits = NULL, col.names
	 = c("lon", "lat"), minsym = "<", label.resolution = 0, labtxt = NULL,
	boxcol = 0, first.color.trans = TRUE, mai = c(0, 1, 0, 1), leftrat = 0.1,
	labbox = TRUE, csi = NULL)
{
     geopar <- getOption("geopar")
     if(!is.null(csi)) cex <- cex*csi/0.12 # Compatibility
	if(!is.null(attributes(grd)$grid)) { 
		z <- grd
		grd <- attributes(grd)$grid
	}
	set <- NA
	fact <- 2
	grd <- Set.grd.and.z(grd, z, NULL, set, col.names)
	# Set data on correct form.
	z <- grd$z
	# Perturb z a little
	z <- z + rnorm(length(z)) * 1e-09
	grd <- grd$grd
	if(is.null(levels)) {
		# changed before cont < 2  
		if(nlevels == 0) nlevels <- 10
		levels <- pretty(range(z, na.rm = TRUE), nlevels)
		levels <- levels[2:(length(levels) - 1)]
	}
	ncont <- length(levels)
	#	Set colors if needed
	if(is.null(col)) {
		if(white) {
			# lowest values white.  
			mincol <- 2
			colors <- c(1:(ncont))
			colors <- floor(mincol + ((colors - 1) * (maxcol - 
				mincol))/(length(colors) - 1))
			colors <- c(0, colors)
		}
		else {
			mincol <- 2
			colors <- c(1:(ncont + 1))
			colors <- floor(mincol + ((colors - 1) * (maxcol - 
				mincol))/(length(colors) - 1))
		}
	}
	else colors <- col
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
	ncont <- nlevels <- length(levels)
	cont <- levels
	# change names of variables
	grd <- extract(grd, z, maxn, limits, col.names = col.names)
	# extract.   
	z <- grd$z
	grd <- grd$grd1
	ind <- c(1:length(z))
	ind <- ind[is.na(z)]
	if(length(ind) > 0) {
		if(fill == 0)
			z[ind] <- -99999
		if(fill == 1)
			z[ind] <- 0
		if(fill == 2)
			z[ind] <- mean(z)
	}
	lon <- grd[[col.names[1]]]
	lat <- grd[[col.names[2]]]
	if(only.positive) {
		# put z<0 to 0
		ind <- c(1:length(z))
		ind <- ind[z < mean(z[z > 0])/1000 & z != -99999]
		z[ind] <- mean(z[z > 0])/1000
	}
	# Check if a setup from geoplot is to be used.  
	cond1 <- col.names[1] == "lon" && col.names[2] == "lat"
	cond2 <- col.names[1] == "x" && col.names[2] == "y" && geopar$
		projection == "none"
	if(cond1 || cond2) {
		oldpar <- selectedpar()
		on.exit(par(oldpar))
		par(geopar$gpar)
		if(geopar$cont)
			par(plt = geopar$contlines)
	}
	if(cex != 0)
		par(cex = cex)
	nx <- length(lon)
	ny <- length(lat)
	mcont <- mean( - cont[1:(ncont - 1)] + cont[2:(ncont)])
	lon1 <- matrix(lon, nx, ny)
	lat1 <- t(matrix(lat, ny, nx))
	# Transform the matrices if lat,lon.  Proj only transforms if col.names=0
	if(col.names[1] == "lon" & col.names[2] == "lat") {
		x1 <- Proj(lat1, lon1, geopar$scale, geopar$b0, geopar$b1,
			geopar$l1, geopar$projection, col.names)
		y1 <- x1$y
		x1 <- x1$x
	}
	else {
		# no projection.
		x1 <- lon1
		y1 <- lat1
	}
	# Check what is inside borders if given.  
	cutreg <- FALSE
	lx <- length(x1)
	x1 <- x1 + rnorm(lx)/1000000
	inni <- 0
	indd <- c(0, 1, 4, 1, 2, 4, 2, 3, 4, 3, 0, 4)
	cont <- c(cont, max(c(max(abs(cont)) * 1.1, max(z) * 1.1)))
	# change.
	cont <- c(min(c(min(z[z != -99999]) - 1, cont[1] - 1)), cont)
	if(cont[2] - cont[1] < 1)
		cont[1] <- cont[2] - 1
	ncont <- ncont + 2
	if(!labels.only) {
		nel <- (nx - 1) * (ny - 1)
		el <- matrix(0, nel, 4)
		err <- 0
		# error
		#    flag
		if(working.space == 0) {
			# Try to calculate working space
			z1 <- z[z > -99998]
			zm <- mean(abs(z1[3:(length(z1) - 1)] - z1[2:(length(
				z1) - 2)]))
			cm <- mean(abs(cont[3:(length(cont) - 1)] - cont[2:
				(length(cont) - 2)]))
			nel <- (nx - 1) * (ny - 1) * 4
			working.space <- round(nel * (zm/cm + 3) * fact * 2)/
				2
			if(working.space <= 20000)
				working.space <- 20000
			# maximum lenght of contour lines
			if(working.space > 200000) working.space <- 10 * length(
					x1)
		}
		polyx <- c(1:working.space)
		polyx[1:working.space] <- 0
		polyy <- c(1:working.space)
		polyy[1:working.space] <- 0
		print(paste("calculated working space is", working.space))
		polyy <- .C("elcont", PACKAGE = "geo", 
			as.double(c(x1)),
			as.double(c(y1)),
			as.double(c(z)),
			as.integer(lx),
			as.integer(nx),
			as.integer(ny),
			as.double(cont),
			as.integer(ncont),
			as.double(polyx),
			as.double(polyy),
			as.integer(el),
			as.integer(nel),
			as.integer(working.space),
			as.integer(err),
			as.integer(inni),
			as.integer(cutreg),
			as.integer(indd),
			as.integer(white))
		err <- polyy[[14]]
		if(err == 1) {
			print("error, working.space not big enough, try calling program again with bigger working.space"
				)
			return(invisible())
		}
		polyx <- polyy[[9]]
		polyy <- polyy[[10]]
		ind <- c(1:length(polyy))
		ind <- ind[polyy > 90000]
		indmax <- ind[length(ind)]
		print(paste(" used working space is ", indmax))
		polyy[ind] <- NA
		col <- polyx[ind]
		# find the color
		polyx[ind] <- NA
		# between polygons.  
		polyx <- polyx[1:indmax]
		polyy <- polyy[1:indmax]
		#		Finding the colors.  
		cont <- cont + 1e-05
		# numerical problem
		ncol <- cut(col, cont,labels=FALSE)  # labels=FALSE R ver.
		ncol <- colors[ncol]
		ind <- c(1:length(ncol))
		ind <- ind[is.na(ncol)]
		ncol[ind] <- 0
		if(bordercheck) {
			eps <- 0.0001
			lx <- geopar$limx - bcrat * (geopar$limx - mean(geopar$
				limx))
			ly <- geopar$limy - bcrat * (geopar$limy - mean(geopar$
				limy))
			lx[1] <- lx[1] - eps
			lx[2] <- lx[2] + eps
			ly[1] <- ly[1] - eps
			ly[2] <- ly[2] + eps
			ind <- c(1:length(polyx))
			ind <- ind[is.na(polyx) | (polyx < lx[2] & polyx > lx[
				1] & polyy < ly[2] & polyy > ly[1])]
			polyx <- polyx[ind]
			polyy <- polyy[ind]
			n <- length(polyx)
			ind <- c(1:length(polyx))
			ind <- ind[!is.na(polyx)]
			if(ind[1] != 1) {
				ind2 <- c(1:(ind[1] - 1))
				polyx <- polyx[ - ind2]
				polyy <- polyy[ - ind2]
				ncol <- ncol[ - ind2]
			}
			ind <- c(1:length(polyx))
			ind <- ind[is.na(polyx)]
			ind1 <- c(1:length(ind))
			ind1 <- ind1[diff(ind) == 1]
			ind2 <- ind[diff(ind) == 1]
			if(length(ind1) > 0) {
				ncol <- ncol[ - (ind1 + 1)]
			}
		}
		polygon(polyx, polyy, col = ncol, border = FALSE)
	}
	# 	Add  labels around plot
	if(length(label.location) == 1) if(label.location == "locator") {
			# use the locator.
			if(cond1 | cond2) label.location <- geolocator(n = 2)
				 else label.location <- locator(n = 2)
		}
	if(length(label.location) > 1) {
		#label located somewhere in drawing
		if(labbox) paint.window(Proj(label.location, col.names = 
				col.names), border = TRUE, col.names = c("y",
				"x"), col = boxcol)
		label.location <- Proj(label.location, scale = geopar$scale,
			b0 = geopar$b0, b1 = geopar$b1, l1 = geopar$l1, 
			projection = geopar$projection, col.names = col.names)
		if(labels == 1) {
			# labels for each contour line.  
			labels1(levels.1, digits, colors.1, xlim = 
				label.location$x, ylim = label.location$y,
				minsym = minsym, label.resolution = 
				label.resolution, labtxt = labtxt, 
				first.color.trans = first.color.trans, mai = 
				mai, leftrat = leftrat)
		}
		else {
			#more of a constant label. 
			labels2(levels.1, digits, colors.1, xlim = 
				label.location$x, ylim = label.location$y)
		}
	}
	if(geopar$cont && labels != 0) {
		# if labels needed.  
		par(plt = geopar$contlab)
		par(new = TRUE)
		plot(c(0, 1, 1, 0, 0), c(0, 0, 1, 1, 0), type = "l", axes = FALSE,
			xlab = " ", ylab = " ")
		if(labels == 1) {
			# labels for each contour line.  
			labels1(levels.1, digits, colors.1, fill = geopar$
				cont, minsym = minsym, label.resolution = 
				label.resolution, labtxt = labtxt, 
				first.color.trans = first.color.trans, mai = 
				mai, leftrat = leftrat)
		}
		else {
			#more of a constant label. 
			labels2(levels.1, digits, colors.1, fill = geopar$
				cont)
		}
	}
	return(invisible())
}

