###########################################################################/**
# @RdocFunction capturePlot
#
# @title "Captures a plot such that it can be redrawn later/elsewhere"
#
# \description{
#   @get "title".
#
#   \emph{This feature is only supported in R (>= 3.3.0).}
# }
#
# @synopsis
#
# \arguments{
#   \item{expr}{The @expression of graphing commands to be evaluated.}
#   \item{envir}{The @environment where \code{expr} should be evaluated.}
#   \item{type}{The type of graphics device used in the background.
#    The choice should not matter since the result should be identical
#    regardless.}
#  \item{...}{Additional arguments passed to the graphics device.}
# }
#
# \value{
#   Returns a \code{recordedplot} object, which can be
#   \code{\link[grDevices]{replayPlot}()}:ed.  If replayed in an
#   interactive session, the plot is displayed in a new window.
#   For conveniency, the object is also replayed when \code{print()}:ed.
# }
#
# \details{
#   Note that plot dimensions/aspect ratios are not recorded.  This
#   means that one does not have to worry about those when recording
#   the plot.  Instead, they are specified when setting up the graphics
#   device(s) in which the recorded plot is replayed (see example).
# }
#
# @examples "../incl/capturePlot.Rex"
#
# @author
#
# \seealso{
#   Internally \code{\link[grDevices]{recordPlot}()} is used.
# }
#
# \references{
#  [1] Paul Murrell et al.,
#      \emph{Recording and Replaying the Graphics Engine Display List},
#      December 2015.
#      \url{https://www.stat.auckland.ac.nz/~paul/Reports/DisplayList/dl-record.html}\cr
# }
#
# @keyword device
#*/###########################################################################
capturePlot <- function(expr, envir=parent.frame(), type=pdf, ...) {
  if (getRversion() < "3.3.0") {
    throw(sprintf("Insufficient R version. R.devices::capturePlot() requires R (>= 3.3.0): ", getRversion()))
  }

  expr <- substitute(expr)

  ## Plot to /dev/null file (or NUL on Windows)
  type(nullfile(), ...)
  on.exit(dev.off())

  dev.control("enable")
  eval(expr, envir=envir)
  recordPlot()
}
