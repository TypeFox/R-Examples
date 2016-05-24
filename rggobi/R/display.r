# View displays
# See methods for details
# 
# @keyword internal 
displays <- function(x) UseMethod("displays", x)

# Get GGobi displays
# Gets list of displays in the specified GGobi instance
# 
# A display basically corresponds to a window in GGobi.  A display
# may contain mutliple plots within it.  For example, the scatterplot
# matrix contains $p * p$ plots.
# 
# Use this function to obtain a reference to a display (they are 
# numbered in the order they are created) so you can change
# display mode, set variables (\code{\link{variables<-.GGobiDisplay}}),
# or save a static image to disk.
# 
# @seealso \code{\link{display}} to create displays 
# @alias displays.GGobi
# @arguments GGobi object
# @keyword dynamic
#X g <- ggobi(mtcars)
#X displays(g)
#X display(g[1])
#X displays(g)
displays.GGobi <- function(x) {
 .GGobiCall("getDisplays", .gobi = x)
}

# Get display dataset
# Returns a link to the GGobiData (dataset) object associated with this display.
# 
# See \code{\link{[.GGobi}} for more information on 
# 
# @arguments GGobiDisplay object
# @arguments ggobi reference
# @keyword manip 
#X g <- ggobi(mtcars)
#X d <- displays(g)[[1]]
#X dataset(d)
dataset.GGobiDisplay <- function(x, .gobi=ggobi(x)) {
 .GGobiCall("getDisplayDataset", x, .gobi=.gobi)
}

# Get display GGobi
# Returns a link to the GGobi object associated with this display
# 
# @arguments GGobiDisplay object
# @keyword internal
ggobi.GGobiDisplay <- function(data, ...) {
	.GGobiCall("getGGobiForDisplay", data)
}

# Print method for GGobiDisplay
# Shows display type (first element in class vector)
# 
# @arguments GGobiDisplay object
# @keyword internal 
print.GGobiDisplay <- function(x, ...) {
	print(class(x)[1])
}

# Close display
# Closes the referenced display.  The R variable will be invalid after this call.
# 
# @arguments GGobiDisplay object to close
# @keyword internal 
#X g <- ggobi(mtcars)
#X displays(g)
#X close(displays(g)[[1]])
#X displays(g)
close.GGobiDisplay <- function(con, ...) {
	.GGobiCall("closeDisplay", con)
}

# Length method for GGobiDisplay
# Returns the number of plots within a given display
# 
# @arguments GGobiDisplay object from which to retrieve the number of plots
# @keyword internal 
length.GGobiDisplay <- function(x) {
	.GGobiCall("getNumPlotsInDisplay",  x)
}


# Save picture of plot (and window) to disk
# This allows you to make a static copy of a GGobiDisplay.
# 
# If you want to make publicaiton quality graphics, you should
# probably use the DescribeDisplay plugin and package.   This
# will recreate a GGobiDisplay in R, and so can produce high-quality
# vector (eg. pdf) output.  See \url{http://www.ggobi.org/describe-display}
# for more information
# 
# @arguments GGobiDisplay to save
# @arguments path to save to
# @arguments type of file to save
# @arguments if TRUE, save only plot, otherwise save surrounding GUI elements as well
# @keyword hplot 
#X g <- ggobi(mtcars)
#X ggobi_display_save_picture(displays(g)[[1]], "test.png")
ggobi_display_save_picture <- function(display=displays(ggobi_get())[[1]], path="ggobi_display.png", filetype="png", plot.only = FALSE) {
	display_widget <- .ggobi_display_get_widget(as.RGtkObject(display))
	if (plot.only) {
          if (inherits(display, "GGobiScatmatDisplay"))
            display_widget <- display_widget[[2]][[1]]
          display_widget <- display_widget[[2]][[3]][["widget"]]
        }
        if (boundGTKVersion() >= "2.14.0")
          drawable <- display_widget$getSnapshot()
        else drawable <- display_widget[["window"]]
        
	widget_size <- drawable$getSize()

	pixbuf <- gdkPixbufGetFromDrawable(src = drawable, cmap = NULL, 
                                           src.x = 0, src.y = 0,
                                           dest.x = 0, dest.y = 0, 
                                           width = widget_size$width,
                                           height = widget_size$height)
	
	pixbuf$save(path, filetype)
        invisible(path)
}

.ggobi_display_get_widget <- function(display) {
  .GGobiCall("getDisplayWidget", display)
}


# =========================================================
# Class methods
# =========================================================

# GGobiDisplay types
# Get list of GGobiDisplay types.  An instance of GGobi must be open
# 
# @arguments GGobi reference
# @seealso \code{\link{ggobi_display_make_type}}
# @keyword internal
ggobi_display_types <- function() {
 .GGobiCall("getDisplayTypes", .gobi = ggobi_get())
}

# Convert plot name to GGobi plot type
# Used to convert between friendly plot name and internal GGobi name.
# 
# @keyword internal 
ggobi_display_make_type <- function(type) {
	if(!inherits(type, "GType")) {
    types <- ggobi_display_types()
    id <- match(type, names(types))
    if(is.na(id)) {
      id <- match(type, sapply(types, pmodes))
    }

    if(is.na(id))
      stop("Unrecognized plot type")

    type <- types[[id]]
  }
	type
}

