#' Draw plot
#'
#' This funtion draws a plot to screen, a file, or both.
#'
#' To send the plot to screen, set \code{device} to NA (default).
#' Optionally, to print the plot on screen to a file, specify \code{file}.
#'
#' If \code{device} is \code{NULL}, the plot will be sent directly to the
#' the specified \code{file} using a printing device inferred from the file 
#' extension (no graphical window will open).
#'
#' Set the global option \code{plot.device} to affect multiple plots.
#' Graphical parameters including \code{width}, \code{height}, \code{res},
#' \code{units} are obtained from the global option \code{getOption("plot")}.
#'
#' @param expr   expression for plotting
#' @param file   filename
#' @param device plot device
#' @param width  plot width [default: 5]
#' @param height plot height [default: 5]
#' @param aspect.ratio  ratio of width to height
#' @param units  unit of plot dimension [default: "in"]
#' @param res    bitmap resolution, used only by bitmap formats [default: 300]
#' @param mkpath   whether to create parent directories
#'                 (if they do not already exists)
#' @param symlink  whether to create a symlink to file with a simplified
#'                 filename (ignored if file is not a \code{filename} object);
#'                 an existing file will not be overwritten but an existing
#'                 symlink will be
#' @param ...    other arguments passed to the plot device function
#' @export
#'
#' @examples
#' \dontrun{
#' # Set device to jpeg (remember to update file extensions for printed plots)
#' options(plot.device=jpeg)
#' qdraw(plot(1:10), "plot.jpeg")
#'
#' # Enable automatic plot format inference
#' options(plot.device=NULL)
#'
#' # Plot directly to file (format is inferred from filename extension)
#' qdraw(plot(1:10), "plot.pdf")
#'
#' # Plot to screen, then print to file (display will not be closed)
#' qdraw(plot(1:10), "plot.png", device=NA)
#' 
#' # If an error occurs, be sure to clear the current plot
#' dev.off()
#' # or clear all plots
#' graphics.off()
#' }
#'
qdraw <- function(
	expr,
	file=NULL, device=getOption("plot.device"),
	width=NULL, height=NULL,
	aspect.ratio=NULL, units=NULL, res=NULL,
	mkpath=TRUE, symlink=TRUE,
	...
) {

	# Set options
	opts <- getOption("plot");
	if (is.null(res)) res <- opts$res;
	if (is.null(units)) units <- opts$units;

	# Set dimensions
	if (is.null(aspect.ratio)) {
		# use default dimensions if not specified
		if (is.null(width)) width <- opts$width;
		if (is.null(height)) height <- opts$height;
	} else {
		# use one dimension and aspect ratio to complet specification
		if (is.null(height)) {
			if (is.null(width)) width <- opts$width;
			height <- width / aspect.ratio;
		} else if (is.null(width)) {
			if (is.null(height)) height <- opts$height;
			width <- height * aspect.ratio;
		}
	}

	filename <- file;

	# Create parent directories to file
	if (!is.null(filename) && mkpath) filenamer::make_path(filename);

	# Convert filenamer::filename to character
	if (is.filename(file)) file <- as.character(file);

	# Infer the device
	if (is.null(device)) {
		# infer device from file extension
		if (!is.null(file)) {
			device <- .infer_device_type(.infer_file_type(file));
		} else {
			# no file is available: plot to screen instead
			device <- NA;
		}
	} else if (is.character(device)) {
		# device is refered to by character string
		device <- .infer_device_type(device);
	}

	# Set up the device
	if (is.function(device)) {

		device.class <- .get_device_class(device);
		if (device.class == "vector") {
			# set oneflie=FALSE to facilitate embedded base plot in grid
			device(width=width, height=height, file=file, onefile=FALSE, ...);
		} else if (device.class == "bitmap") {
			device(width=width, height=height, filename=file,
				units=units, res=res, ...);
		} else {
			# Assume device plots to screen (no printing later)
			device(width=width, height=height, ...);
		}
		# Set file to NULL to prevent dev.print
		file <- NULL;

	} else if (is.na(device)) {

		# Use default device (i.e. screen)
		dev.new(width=width, height=height, ...);

	} else {
		stop("Unknown plot output device type")
	}

	# Evaluate the plotting expression
	g <- eval(substitute(expr), parent.frame(), parent.frame());
	if (inherits(g, "ggplot")) print(g);

	# Finalize the plotting device
	if (is.null(device) || is.function(device) || is.character(device)) {
		# Assume plot was sent to a file: close device
		# NB  Do not set plot.device option to a screen device!
		dev.off();
	} else if (is.na(device)) {
		# plot was sent to screen
		if (!is.null(file)) {
			# plot.print is specified: print to file
			# infer file and device type from file extension
			print.device <- .infer_device_type(.infer_file_type(file));
			device.class <- .get_device_class(print.device);
			if (device.class == "vector") {
				# set oneflie=FALSE to facilitate embedded base plot in grid
				dev.print(print.device, width=width, height=height, file=file,
					onefile=FALSE, ...);
			} else if (device.class == "bitmap") {
				dev.print(print.device, width=width, height=height, filename=file,
					units=units, res=res, ...);
			} else {
				stop("Could not infer plot output device type")
			}
		}
		# otherwise, leave device open and do nothing
	} else {
		stop("Unknown plot output device type")
	}

	# Create symbolic link
	if (!is.null(filename) && symlink) .symlink_simplified(filename);

	invisible()
}

.get_device_class <- function(device) {
	if (
		# vector formats
		identical(device, pdf) ||
		identical(device, postscript) ||
		identical(device, cairo_pdf) ||
		identical(device, cairo_ps) ||
		identical(device, svg)
	) {
		"vector"

	} else if (
		# bitmap formats
		identical(device, png) ||
		identical(device, jpeg) ||
		identical(device, bmp) ||
		identical(device, tiff)
	) {
		"bitmap"

	} else {
		NA
	}
}

.infer_device_type <- function(name) {
	# standardize type names, as necessary
	name <- tolower(name);
	device.type <- switch(name,
		jpg = "jpeg",
		ps = "postscript",
		name
	);
	if (device.type %in%
		c(
			"pdf", "postscript", "cairo_pdf", "cairo_ps", "svg",
			"png", "jpeg", "bmp", "tiff"
		)
	) {
		# device type one of the known devices
		return( get(device.type) )
	} else {
		stop("Could not infer plot output device type")
	}
}
