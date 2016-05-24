Cairo <- function(width = 7, height = 7, pointsize = 12,
                  surface = c("screen", "png", "pdf", "ps", "svg"),
                  filename = NULL)
{
  if (missing(surface) && !missing(filename))
    surface <- tolower(gsub(".*\\.", "", basename(filename)))
  surface <- match.arg(surface)
  if (missing(filename) && surface != "screen")
    filename <- paste("Rplots", surface, sep=".")
  surface_info <- c(surface, filename)
  .C("do_Cairo", as.numeric(width), as.numeric(height), as.numeric(pointsize),
    as.character(surface_info), PACKAGE="cairoDevice")
  return(invisible(TRUE))
}

Cairo_png <- function(filename, width = 7, height = 7, pointsize = 12)
{
  Cairo(width, height, pointsize, "png", filename)
}
Cairo_pdf <- function(filename, width = 7, height = 7, pointsize = 12)
{
  Cairo(width, height, pointsize, "pdf", filename)
}
Cairo_ps <- function(filename, width = 7, height = 7, pointsize = 12)
{
  Cairo(width, height, pointsize, "ps", filename)
}
Cairo_svg <- function(filename, width = 7, height = 7, pointsize = 12)
{
  Cairo(width, height, pointsize, "svg", filename)
}

asCairoDevice <- function(widget, pointsize = 12, width = 500, height = 500)
{
  w <- -1
  h <- -1
  if (inherits(widget, "GtkPrintContext")) {
    w <- widget$getWidth()
    h <- widget$getHeight()
    widget <- widget$getCairoContext()
  } else if (inherits(widget, "Cairo")) {
    stopifnot(width >= 0)
    stopifnot(height >= 0)
    w <- as.numeric(width)
    h <- as.numeric(height)
  } else if (!inherits(widget, "GtkDrawingArea") &&
             !inherits(widget, "GdkDrawable"))
    {
      stop("Object being used as a Cairo Device must be a GtkDrawingArea, ",
           "GdkDrawable, GtkPrintContext or Cairo context")
    }
  if (!inherits(widget, "Cairo") && (!missing(width) || !missing(height)))
    stop("'width' and 'height' ignored unless 'widget' is a Cairo context")
  .Call("do_asCairoDevice", widget, as.numeric(pointsize), w, h,
        PACKAGE="cairoDevice")
}
