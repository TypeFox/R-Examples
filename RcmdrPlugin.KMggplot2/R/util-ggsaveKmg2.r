#' A Modified Function of \code{ggsave} for \code{recordedplot} Class Objects
#'
#' This is a modified function of \code{ggsave} for \code{recordedplot} class
#' objects.
#'
#' @param filename file name/filename of plot
#' @param plot plot to save, defaults to last plot displayed
#' @param device device to use, automatically extract from file name extension
#' @param path path to save plot to (if you just want to set path and not filename)
#' @param scale scaling factor
#' @param width width (defaults to the width of current plotting window)
#' @param height height (defaults to the height of current plotting window)
#' @param units units for width and height when either one is explicitly specified (in, cm, or mm)
#' @param dpi dpi to use for raster graphics
#' @param ... other arguments passed to graphics device
#' @seealso \code{\link[ggplot2:ggsave]{ggsave}}
#'
#' @rdname util-ggsaveKmg2
#' @keywords hplot
#' @importFrom ggplot2 last_plot
#' @importFrom grDevices dev.off
#' @importFrom graphics par
#' @importFrom utils capture.output
#' @export
ggsaveKmg2 <- function(
  filename = default_name(plot), plot = last_plot(), device = default_device(filename),
  path = NULL, scale = 1, width = par("din")[1], height = par("din")[2], units = c("in", "cm", "mm"),
  dpi = 300, ...) {

  if (!inherits(plot, "ggplot") && !inherits(plot, "recordedplot")) 
    stop("plot should be a ggplot2 plot or a recordedplot plot")

  eps <- ps <- function(..., width, height)  
    grDevices::postscript(..., width = width, height = height, onefile = FALSE,
      horizontal = FALSE, paper = "special")
  tex <- function(..., width, height) 
    grDevices::pictex(..., width = width, height = height)
  pdf <- function(..., version = "1.4") 
    grDevices::pdf(..., version = version)
  svg <- function(...) 
    grDevices::svg(...)
  wmf <- function(..., width, height) 
    grDevices::win.metafile(..., width = width, height = height)

  png <- function(..., width, height) 
    grDevices::png(...,  width = width, height = height, res = dpi, units = "in")
  jpg <- jpeg <- function(..., width, height) 
    grDevices::jpeg(..., width = width, height = height, res = dpi, units = "in")
  bmp <- function(..., width, height) 
    grDevices::bmp(...,  width = width, height = height, res = dpi, units = "in")
  tiff <- function(..., width, height) 
    grDevices::tiff(..., width = width, height = height, res = dpi, units = "in")
  
  default_name <- function(plot) { 
    paste("default_plot.pdf", sep = "")
  }
  
  default_device <- function(filename) {
    pieces <- strsplit(filename, "\\.")[[1]]
    ext <- tolower(pieces[length(pieces)])
    match.fun(ext)
  }

  units <- match.arg(units)
  convert_to_inches <- function(x, units) {
    x <- switch(units,
      `in` = x,
      cm = x / 2.54,
      mm = x / 2.54 /10
    )
  }

  convert_from_inches <- function(x, units) {
    x <- switch(units,
      `in` = x,
      cm = x * 2.54,
      mm = x * 2.54 * 10
    )
  }

  # dimensions need to be in inches for all graphic devices
  # convert width and height into inches when they are specified
  if (!missing(width)) {
    width <- convert_to_inches(width, units)
  }
  if (!missing(height)) {
    height <- convert_to_inches(height, units)
  }
  # if either width or height is not specified, display an information message
  # units are those specified by the user
  if (missing(width) || missing(height)) {
    message("Saving ", prettyNum(convert_from_inches(width * scale, units), digits = 3),
            " x ", prettyNum(convert_from_inches(height * scale, units), digits = 3),
            " ", units, " image")
  }

  width <- width * scale
  height <- height * scale
  
  if (!is.null(path)) {
    filename <- file.path(path, filename)
  }
  device(file = filename, width = width, height = height, ...)
  on.exit(capture.output(dev.off()))
  print(plot)
  
  invisible()

}
