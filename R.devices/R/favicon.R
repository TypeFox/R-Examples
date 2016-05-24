#########################################################################/**
# @RdocFunction favicon
#
# @title "Favicon graphics device"
#
# \description{
#  Device driver for PNG favicons (Favorite icon) also known as
#  shortcut icon, Web site icon, tab icon or bookmark icon.
#  This driver is the same as the png driver where some arguments
#  have different default values.
# }
#
# @synopsis
#
# \arguments{
#   \item{file}{Default file name (pattern).}
#   \item{width, height}{The width and height of the figure.}
#   \item{par}{A named @list of graphical parameters to use.}
#   \item{...}{Other arguments accepted by \code{png()}.}
# }
#
# \value{
#   A plot device is opened; nothing is returned.
# }
#
# \examples{\dontrun{
#   favicon(width=32L)
#
#   # is identical to
#
#   suppressWarnings({
#     png("favicon.png", width=32L, height=32L, bg="transparent",
#                        par=list(mar=c(0,0,0,0)))
#   })
# }}
#
# \keyword{device}
#
# \seealso{
#   Internally, @see "grDevices::png" is used.
#   %% Add HTML Favicon script
#   %% toFavicon(plot(1, col="blue", bg="yellow",
#   %%           pch=21, cex=4, lwd=4, axes=FALSE))
#   \if{html}{\out{
#    <script>
#     var link = document.createElement('link');
#     link.rel = 'icon';
#     link.href = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAMAAABEpIrGAAAADFBMVEX9/v0AAP/9/v3//wBEQjoBAAAABHRSTlP//wD//gy7CwAAAGJJREFUOI3N0rESwCAIA9Ag///PXdoiBk0HhmbNO49DMETQCexNCSyFgdlGoO5DYOr9ThLgPosA7osIQP0sHuDOog8UI/ALa988wzdwXJRctf4s+d36YPTJ6aMd8ux3+QO4ABTtB85yDAh9AAAAAElFTkSuQmCC";
#     document.getElementsByTagName('head')[0].appendChild(link);
#    </script>
#   }}
# }
#
# @author
#
# @keyword device
# @keyword internal
#*/#########################################################################
favicon <- function(filename="favicon.png", width=32L, height=width, bg="transparent", par=list(mar=c(0,0,0,0)), ...) {
  # Argument 'width' and 'height':
  if (height != width) {
    throw("The width and the height must be the same for a favicon: ",
          width, " != ", height)
  }

  # Argument 'par':
  if (!is.null(par)) {
    if (!is.list(par) || is.null(names(par))) {
      throw("Argument 'par' has to be a named list: ", mode(par))
    }
  }

  # Create PNG file
  # png() will generate a warning that "width=.. , height=.., are unlikely
  # values in pixels" if width < 20.
  suppressWarnings({
    png(filename=filename, width=width, height=height, bg=bg, ...)
  })

  # Set graphical parameters
  par(par)
} # favicon()


############################################################################
# HISTORY:
# 2014-09-15
# o Created.
############################################################################

