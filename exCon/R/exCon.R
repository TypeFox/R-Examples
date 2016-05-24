##' Explore Contour Data Interactively
##'
##' These functions compute contour lines from matrix data and display them
##' in an interactive web page using the d3 javascript library. \code{exCon}
##' displays slices along the x and y directions.  \code{exCon2} does not
##' display the slices and is faster.
##' 
##' @param M A matrix.
##'
##' @param x A vector of numeric values giving the locations of the grid defining
##' the matrix.  Must have length \code{nrow(M)}.  These values become the scale
##' along the x axis.
##'
##' @param y A vector of numeric values giving the locations of the grid defining
##' the matrix.  Must have length \code{ncol(M)}. These values become the scale
##' along the y axis.
##'
##' @param nlevels  Integer.  The number of contour levels desired.  Ignored if \code{levels}
##' is given.
##'
##' @param levels  Numeric.  A vector of values (altitudes if you will) at which to
##' compute the contours.
##'
##' @param browser Character.  Something that will make sense to your OS.  Only
##' necessary if you want to overide your system specified browser as understood by
##' \code{R}.  See below for further details.
##'
##' @param minify Logical.  Shall the JavaScript be minified?  This improves
##' performance.  However, it requires package \code{js} which in turn requires
##' package \code{V8}.  The latter is not available on all platforms.  Details
##' may be available at \url{https://github.com/jeroenooms/v8}
##'
##' @return The path to the temporary directory containing the web page files.
##' The side effect is an interactive web page.  The temporary directory
##' is deleted when you quit R, but you can use the return value to
##' save the files in a different location.
##'
##' @section Details: The computation of the contour lines is handled by
##' \code{\link[grDevices]{contourLines}}.  The result here, however, is transposed so that the
##' output has the same orientation as the original matrix. This is necessary because
##' \code{\link[graphics]{contour}} tranposes its output: "Notice that
##' \code{contour} interprets the \code{z} matrix as a table of 
##' \code{f(x[i], y[j])} values, so that the x axis corresponds to row number
##' and the y axis to column number, with column 1 at the bottom, i.e. a 90 degree
##' counter-clockwise rotation of the conventional textual layout."
##' 
##' @section Interpretation:  The contour lines are an interpolation of the data
##' in the matrix.  In \code{exCon}, the slices are the actual values in the matrix
##' row or column
##' connected point-to-point.  Thus a maximum in a slice may not correspond to 
##' a peak in the contour plot.
##' 
##' @section Browser Choice: The browser is called by
##' \code{\link[utils]{browseURL}}, which
##' in turn uses \code{options("browser")}.  Exactly how this is handled
##' is OS dependent.
##'
##' @section RStudio Viewer: If browser is \code{NULL}, you are using RStudio, and a viewer is specified, this will be called.  You can stop this by with \code{options(viewer = NULL)}.
##'
##'
##' @section Browser Choice (Mac): On a Mac, the default browser is called
##' by \code{/bin/sh/open}
##' which in turn looks at which browser you have set in the system settings.  You can
##' override your default with
##' \code{browser = "/usr/bin/open -a 'Google Chrome'"} for example.
##' Testing shows that on a Mac, Safari and Chrome perform correctly,
##' but in Firefox the mouse cursor is slightly offset from the guides.  However,
##' the slices chosen and the values displayed are correct.
##'
##' @section Browser Choice (Other Systems):  \code{exCon} has been tested
##' on a Windows 7
##' professional instance running in VirtualBox using Firefox and Chrome, and
##' runs correctly (Firefox has the same mouse position issue as mentioned above).
##' It runs similarly on Window 8.
##'
##' @section Browser Choice & Performance:  You can check the performance of
##' your browser at peacekeeper.futuremark.com  The most relevant score for
##' exCon is the rendering category.  In limited testing, Chrome does the best.
##'
##' @section Performance Limits (YMMV): On a 4-year old MacBook Pro, with 8 Gb
##' RAM and an Intel Core i7 chip, a
##' 1500 x 1500 matrix with 1 contour level 
##' requires about 30 seconds for R to render the web page using Chrome, Safari
##' or Firefox. A 2K x 2K matrix appears to be too large to handle.
##' Approximately the same result is obtained using Windows 8 and Chrome.
##'
##' @name exCon
##' @aliases exCon2 exCon
##' @rdname exCon
##' @export exCon exCon2
##' @import jsonlite
##' @importFrom utils browseURL
##' @importFrom grDevices contourLines
##' @keywords plot interactive
##'
##' @examples
##' if (interactive()) {
##' require(jsonlite)
##'
##' # minify is FALSE in the examples as not all platforms support the required pkgs (see above)
##'
##' exCon(M = volcano, minify = FALSE)
##'
##' exCon2(M = volcano, minify = FALSE) # no slices
##'
##' # This next example will label the axes with the actual values, relative to the
##' # lower left corner (original data collected on 10 meter grid).  Giving
##' # x and y affects only the scale, and the native values displayed at the top.
##' 
##'  exCon(M = volcano, minify = FALSE,
##'  x = seq(from = 0, by = 10, length.out = nrow(volcano)),
##'  y = seq(from = 0, by = 10, length.out = ncol(volcano)))
##' }
##'
exCon <- function(M = NULL,
	x = seq(0, 1, length.out = nrow(M)),
	y = seq(0, 1, length.out = ncol(M)),
	nlevels = 5,
	levels = pretty(range(M, na.rm = TRUE), nlevels),
	browser = NULL, minify = TRUE) {

	# Bryan A. Hanson, DePauw University, April 2014
	# This is the R front end controlling everything
	# The html and js files called indirectly written by
	# Kristina Mulry and Bryan A. Hanson

	if (is.null(M)) stop("You must provide a matrix M")
	if (!is.matrix(M)) stop("M must be a matrix")

	# M is the raw, topographic data matrix:
	# think of it as altitudes on an x,y grid

	# Major steps
	# 1. Convert M to contour lines (CL)
	# 2. Read in the existing js, then pre-pend
	#    M and CL to to it and write it back out
	# 3. Call a browser on the html to open it automatically

	# Compute contour lines
	# (Eventually need to supply nlevels & levels,
	# Also a name for the data set)
	dimnames(M) <- list(rep("x", nrow(M)), rep("y", ncol(M)))
	CL <- contourLines(z = M, nlevels = nlevels, levels = levels,
		x = x, y = y)

	# Get the contour lines into JSON format

	CL <- toJSON(CL)

	# Get M into JSON format
	# [ [row1...], [row2...], [lastRow...] ]

	M <- toJSON(t(M))

	# We need the first and last values of x and y
	# to be the xD and yD (the domains)
	# These are in the units of whatever is supplied (native)

	DX <- toJSON(c(x[1], x[length(x)]))
	DY <- toJSON(c(y[1], y[length(y)]))

	# Prepare the data
	
	data1 <- paste("var CL = ", CL, sep = "")
	data2 <- paste("var M = ", M, sep = "")
	data3 <- paste("var Dx = ", DX, sep = "")
	data4 <- paste("var Dy = ", DY, sep = "")

	# Get the JavaScript modules & related files
	
	td <- tempdir()
	fd <- system.file("extdata", package = "exCon")
	eCfiles <- c("eC.css", "eC_globals.js", "eC_controls.js", "eC_contours.js",
	"eC_brushNguides.js", "eC_slices.js", "eC_main.js", "exCon.html")	
	chk2 <- file.copy(from=file.path(fd, eCfiles), to=file.path(td, eCfiles),
		overwrite = TRUE)
	if (!all(chk2)) stop("Copying to temporary directory failed")

	js1 <- readLines(con = file.path(td,"eC_globals.js"))
	js2 <- readLines(con = file.path(td,"eC_controls.js"))
	js3 <- readLines(con = file.path(td,"eC_contours.js"))
	js4 <- readLines(con = file.path(td,"eC_brushNguides.js"))
	js5 <- readLines(con = file.path(td,"eC_slices.js"))
	js6 <- readLines(con = file.path(td,"eC_main.js"))

	scopeFunHeader <- "(function() {"
	scopeFunTail <- "})();"

	# Now write

	# text = c(scopeFunHeader, data1, data2, data3, data4,
		# js1, js2, js3, js4, js5, js6, scopeFunTail)

	text = c(data1, data2, data3, data4,
		js1, js2, js3, js4, js5, js6)

	if (minify) {
		if (requireNamespace("js", quietly = TRUE)) {
			text <- js::uglify_optimize(text, unused = FALSE)
			}
		if (!requireNamespace("js", quietly = TRUE)) {
			stop("You need to install package js to minify the JavaScript code.  See ?exCon for details.")
			}
		}
	
	writeLines(text, sep  = "\n", con = file.path(td,"exCon.js"))

	# Open the file in a browser

	pg <-  file.path(td,"exCon.html")
	if (!is.null(browser)) {
	    browseURL(pg, browser = browser)
		} else {
		# open in RStudio if viewer is not null
	    # similar to htmltools::html_print
			viewer <- getOption("viewer")
		  	if (is.null(browser) && !is.null(viewer)) {
	      		viewer(pg)
		  		} else {
		    		browseURL(pg)
		  			}
		}
  
	# message("The exCon web page is in the following temp directory which is deleted when you quit R: ")
	# message(td)
	invisible(td)
}



exCon2 <- function(M = NULL,
	x = seq(0, 1, length.out = nrow(M)),
	y = seq(0, 1, length.out = ncol(M)),
	nlevels = 5,
	levels = pretty(range(M, na.rm = TRUE), nlevels),
	browser = NULL, minify = TRUE) {

	# Bryan A. Hanson, DePauw University, June 2014
	# This is the R front end controlling everything
	# The html and js files called indirectly written by
	# Kristina Mulry and Bryan A. Hanson

	if (is.null(M)) stop("You must provide a matrix M")
	if (!is.matrix(M)) stop("M must be a matrix")

	# M is the raw, topographic data matrix:
	# think of it as altitudes on an x,y grid

	# Major steps
	# 1. Convert M to contour lines (CL)
	# 2. Read in the existing js, then pre-pend
	#    M and CL to to it and write it back out
	# 3. Call a browser on the html to open it automatically

	# Compute contour lines
	# (Eventually need to supply nlevels & levels,
	# Also a name for the data set)
	dimnames(M) <- list(rep("x", nrow(M)), rep("y", ncol(M)))
	CL <- contourLines(z = M, nlevels = nlevels, levels = levels,
		x = x, y = y)

	# Get the contour lines into JSON format

	CL <- toJSON(CL)

	# Get M into JSON format
	# [ [row1...], [row2...], [lastRow...] ]

	M <- toJSON(t(M))

	# We need the first and last values of x and y
	# to be the xD and yD (the domains)
	# These are in the units of whatever is supplied (native)

	DX <- toJSON(c(x[1], x[length(x)]))
	DY <- toJSON(c(y[1], y[length(y)]))

	# Prepare the data
	
	data1 <- paste("var CL = ", CL, sep = "")
	data2 <- paste("var M = ", M, sep = "")
	data3 <- paste("var Dx = ", DX, sep = "")
	data4 <- paste("var Dy = ", DY, sep = "")

	# Get the JavaScript modules & related files
	
	td <- tempdir()
	fd <- system.file("extdata", package = "exCon")
	eCfiles <- c("eC.css", "eC2_globals.js", "eC2_contours.js",
	"eC2_brushNguides.js", "eC2_main.js", "exCon2.html")	
	chk2 <- file.copy(from=file.path(fd, eCfiles), to=file.path(td, eCfiles),
		overwrite = TRUE)
	if (!all(chk2)) stop("Copying to temporary directory failed")

	js1 <- readLines(con = file.path(td,"eC2_globals.js"))
	js2 <- readLines(con = file.path(td,"eC2_contours.js"))
	js3 <- readLines(con = file.path(td,"eC2_brushNguides.js"))
	js4 <- readLines(con = file.path(td,"eC2_main.js"))

	scopeFunHeader <- "(function() {"
	scopeFunTail <- "})();"

	# Now write

	text = c(scopeFunHeader, data1, data2, data3, data4,
		js1, js2, js3, js4, scopeFunTail)

	# text = c(data1, data2, data3, data4,
		# js1, js2, js3, js4)

	if (minify) {
		if (requireNamespace("js", quietly = TRUE)) {
			text <- js::uglify_optimize(text, unused = FALSE)
			}
		if (!requireNamespace("js", quietly = TRUE)) {
			stop("You need to install package js to minify the JavaScript code.  See ?exCon for details.")
			}
		}
	
	writeLines(text, sep  = "\n", con = file.path(td,"exCon2.js"))

	# Open the file in a browser

	pg <-  file.path(td,"exCon2.html")
	if (!is.null(browser)) {
	    browseURL(pg, browser = browser)
		} else {
		# open in RStudio if viewer is not null
	    # similar to htmltools::html_print
			viewer <- getOption("viewer")
		  	if (is.null(browser) && !is.null(viewer)) {
	      		viewer(pg)
		  		} else {
		    		browseURL(pg)
		  			}
		}
  
	# message("The exCon2 web page is in the following temp directory which is deleted when you quit R: ")
	# message(td)
	invisible(td)
}
