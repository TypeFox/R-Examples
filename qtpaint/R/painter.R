## R wrappers around the Qanviz::Painter API
## These are the lowest level wrappers

## should be called 'qtransform', but qtbase already uses that for a
## convenience constructor (would rather not make it generic)

##' Access the user to device coordinate transform of the paint context.
##'
##' @title User to device transform
##' @param x The paint context
##' @param value The desired user to device transform, as a C++
##' \code{QTransform} object, or a logical scalar for
##' \code{qdeviceTransformEnabled}.
##' @return The user to device transform, as a C++ \code{QTransform} object
##' @author Michael Lawrence
##' @rdname transform
qdeviceTransform <- function(x) {
  stopifnot(inherits(x, "Painter"))
  .Call("qt_qtransform_Painter", x)
}
##' @rdname transform
##' @usage qdeviceTransform(x) <- value
`qdeviceTransform<-` <- function(x, value) {
  stopifnot(inherits(x, "Painter"))
  stopifnot(inherits(value, "QTransform"))
  .Call("qt_qsetTransform_Painter", x, value, PACKAGE = "qtpaint")
}

##' @rdname transform
##' @usage qdeviceTransformEnabled(x) <- value
`qdeviceTransformEnabled<-` <- function(x, value) {
  stopifnot(inherits(x, "Painter"))
  invisible(.Call("qt_qsetTransformEnabled_Painter", x, as.logical(value)))
}

##' Functions for controlling the stroke and fill colors, including
##' the ability to disable or enable stroking and filling.
##'
##' @title Stroke and fill colors
##' @param x The paint context
##' @param value The stroke or fill color, or a logical value for
##' \code{qHasStroke<-} and \code{qHasFill<-}. A stroke and fill color
##' should either be a C++ \code{QColor} object, a matrix like that
##' returned by \code{\link{col2rgb}} or something coercible to one,
##' like a color name.
##' @return The stroke or fill color, or a logical value for
##' \code{qHasStroke} and \code{qHasFill}.
##' @rdname stroke-fill
##' @author Michael Lawrence
##' @usage qhasStroke(x) <- value
`qhasStroke<-` <- function(x, value) {
  stopifnot(inherits(x, "Painter"))
  invisible(.Call("qt_qsetHasStroke_Painter", x, as.logical(value)))
}

##' @rdname stroke-fill
##' @usage qhasFill(x) <- value
`qhasFill<-` <- function(x, value) {
  stopifnot(inherits(x, "Painter"))
  invisible(.Call("qt_qsetHasFill_Painter", x, as.logical(value)))
}

.normArgStroke <- function(p, color, len) {
  if (length(color) == 1) {
    qstrokeColor(p) <- color
    NULL
  } else .normArgColor(color, len)
}

.normArgFill <- function(p, color, len) {
  if (length(color) == 1) {
    qfillColor(p) <- color
    NULL
  } else .normArgColor(color, len)
}

.normArgColor <- function(color, len) {
  if (is.null(color))
    return(NULL)
  if (inherits(color, "QColor"))
    color <- as.matrix(color) # simplifies C code
  if (!is.matrix(color) || !is.integer(color) || nrow(color) != 4)
    color <- col2rgb(color, TRUE)
  if (!missing(len)) # might drop to vector here, much faster, C code is OK
    color <- recycleVector(color, 4L*len)
  color
}

##' @rdname stroke-fill
##' @usage qstrokeColor(x) <- value
`qstrokeColor<-` <- function(x, value) {
  stopifnot(inherits(x, "Painter"))
  color <- .normArgColor(value)
  if (is.null(color))
    stop("Cannot set a 'NULL' stroke color, use 'NA' for transparent")
  invisible(.Call("qt_qsetStrokeColor_Painter", x, color))
}

##' @rdname stroke-fill
##' @usage qfileColor(x) <- value
`qfillColor<-` <- function(x, value) {
  stopifnot(inherits(x, "Painter"))
  color <- .normArgColor(value)
  if (is.null(color))
    stop("Cannot set a 'NULL' fill color, use 'NA' for transparent")
  invisible(.Call("qt_qsetFillColor_Painter", x, color))
}

##' Sets the font.
##'
##' @title Fonts
##' @param x The paint context
##' @param value The font, as a C++ \code{QFont} object
##' @author Michael Lawrence
##' @name fontReplace
##' @usage qfont(x) <- value
`qfont<-` <- function(x, value)
{
  stopifnot(inherits(x, "Painter"))
  stopifnot(inherits(value, "QFont"))
  invisible(.Call("qt_qsetFont_Painter", x, value))
}

##' Sets the line width
##'
##' @title Line width
##' @param x The paint context
##' @param value A numeric value indicating the line width, in device coords
##' @author Michael Lawrence
##' @name lineWidthReplace
##' @usage qlineWidth(x) <- value
`qlineWidth<-` <- function(x, value) {
  stopifnot(inherits(x, "Painter"))
  invisible(.Call("qt_qsetLineWidth_Painter", x, as.integer(value)))
}

##' Sets the dash pattern
##'
##' @title Dash patterns
##' @param x The paint context
##' @param value A numeric vector indicating the pattern; each element
##' is the length of the corresponding segment
##' @name dashReplace
##' @author Michael Lawrence
##' @usage qdash(x) <- value
`qdash<-` <- function(x, value) {
  stopifnot(inherits(x, "Painter"))
  invisible(.Call("qt_qsetDashes_Painter", x, as.numeric(value)))
}

##' Sets the glyph expansion, equivalent to \code{cex} in base R
##'
##' @title Glyph expansion
##' @param x The paint context
##' @param value Floating point multiplier of the glyph size
##' @author Michael Lawrence
##' @usage qglyphExpansion(x) <- value
`qglyphExpansion<-` <- function(x, value) {
  stopifnot(inherits(x, "Painter"))
  invisible(.Call("qt_qsetGlyphExpansion_Painter", x, as.numeric(value)))
}

##' Enables or disables antialiasing
##'
##' @title Antialiasing
##' @param x The paint context
##' @param value A logical indicating whether antialiasing is enabled
##' @author Michael Lawrence
##' @name antialiasReplace
##' @usage qantialias(x) <- value
`qantialias<-` <- function(x, value) {
  stopifnot(inherits(x, "Painter"))
  invisible(.Call("qt_qsetAntialias_Painter", x, as.logical(value)))
}

##' These functions constitute the primary drawing API. There is
##' support for drawing points, polylines, segments, circles, rectangles,
##' polygons, vector paths, text, images and plot glyphs.
##'
##' @title Drawing API
##' @param p The paint context
##' @param x The X coordinate vector, recycled. For polygons and
##' polylines, NA values separate the graphical primitives.
##' @param y The Y coordinate vector, recycled. For polygons and
##' polylines, NA values separate the graphical primitives.
##' @param stroke The vector of stroke colors, either a C++
##' \code{QColor} object, a matrix returned by \code{\link{col2rgb}}
##' or any valid input to \code{col2rgb}, recycled, or \code{NULL} to
##' disable stroking. Recycled to match the number of primitives.
##' @author Michael Lawrence
##' @rdname painting
qdrawLine <- function(p, x, y, stroke = NULL) {
  stopifnot(inherits(p, "Painter"))
  m <- max(length(x), length(y))
  x <- recycleVector(x, m)
  y <- recycleVector(y, m)
  n <- sum(is.na(x)) + 1L
  invisible(.Call("qt_qdrawPolyline_Painter", p, as.numeric(x), as.numeric(y),
                  .normArgStroke(p, stroke, n)))
}

##' @param x0 The vector of first X coordinates, recycled
##' @param y0 The vector of first Y coordinates, recycled
##' @param x1 The vector of second X coordinates, recycled
##' @param y1 The vector of second Y coordinates, recycled
##' @rdname painting
qdrawSegment <- function(p, x0, y0, x1, y1, stroke = NULL) {
  stopifnot(inherits(p, "Painter"))
  m <- max(length(x0), length(y0), length(x1), length(y1))
  x0 <- recycleVector(x0, m)
  y0 <- recycleVector(y0, m)
  x1 <- recycleVector(x1, m)
  y1 <- recycleVector(y1, m)
  invisible(.Call("qt_qdrawSegments_Painter", p, as.numeric(x0), as.numeric(y0),
                  as.numeric(x1), as.numeric(y1),
                  .normArgStroke(p, stroke, m)))
}

##' @rdname painting 
qdrawPoint <- function(p, x, y, stroke = NULL) {
  stopifnot(inherits(p, "Painter"))
  m <- max(length(x), length(y))
  x <- recycleVector(x, m)
  y <- recycleVector(y, m)
  invisible(.Call("qt_qdrawPoints_Painter", p, as.numeric(x), as.numeric(y),
                  .normArgStroke(p, stroke, m)))
}

##' @param fill The vector of fill colors, either a C++ \code{QColor}
##' object, a matrix returned by \code{\link{col2rgb}} or any valid
##' input to \code{col2rgb}, recycled, or \code{NULL} to disable
##' filling. Recycled to match the number of primitives.
##' @param xleft The vector of left X coordinates for a rectangle, recycled
##' @param ybottom The vector of bottom Y coordinates for a rectangle, recycled
##' @param xright The vector of right X coordinates for a rectangle, recycled
##' @param ytop The vector of top Y coordinates for a rectangle, recycled
##' @rdname painting
qdrawRect <- function(p, xleft, ybottom, xright, ytop, stroke = NULL,
                      fill = NULL)
{
  stopifnot(inherits(p, "Painter"))
  m <- max(length(xleft), length(ybottom), length(xright), length(ytop))
  xleft <- recycleVector(xleft, m)
  ybottom <- recycleVector(ybottom, m)
  xright <- recycleVector(xright, m)
  ytop <- recycleVector(ytop, m)
  invisible(.Call("qt_qdrawRectangles_Painter", p, as.numeric(xleft),
                  as.numeric(ybottom), as.numeric(xright - xleft),
                  as.numeric(ytop - ybottom), .normArgStroke(p, stroke, m),
                  .normArgFill(p, fill, m)))
}

##' @param r The radius of the circle, in device coordinates, recycled
##' @rdname painting 
qdrawCircle <- function(p, x, y, r, stroke = NULL, fill = NULL) {
  stopifnot(inherits(p, "Painter"))
  m <- max(length(x), length(y), length(r))
  x <- recycleVector(x, m)
  y <- recycleVector(y, m)
  r <- recycleVector(r, m)
  invisible(.Call("qt_qdrawCircle_Painter", p, as.numeric(x), as.numeric(y),
                  as.integer(r), .normArgStroke(p, stroke, m),
                  .normArgFill(p, fill, m)))
}

##' @rdname painting 
qdrawPolygon <- function(p, x, y, stroke = NULL, fill = NULL) {
  stopifnot(inherits(p, "Painter"))
  m <- max(length(x), length(y))
  x <- recycleVector(x, m)
  y <- recycleVector(y, m)
  n <- sum(is.na(x)) + 1L
  invisible(.Call("qt_qdrawPolygon_Painter", p, as.numeric(x), as.numeric(y),
                  .normArgStroke(p, stroke, n), .normArgFill(p, fill, n)))
}

##' @rdname painting 
##' @param path A C++ \code{QPainterPath} object describing the glyph,
##' or a list of such objects for \code{qdrawPath}.
qdrawPath <- function(p, path, stroke = NULL, fill = NULL) {
  stopifnot(inherits(p, "Painter"))
  if (inherits(path, "QPainterPath"))
    path <- list(path)
  else path <- as.list(path)
  m <- length(path)
  invisible(.Call("qt_qdrawPath_Painter", p, path,
                  .normArgStroke(p, stroke, m), .normArgFill(p, fill, m)))
}

## Text drawing: a mess

## It seems that (at least in base R graphics) there are three
## different ways to align text: left/bottom, center, right/top.
## For horizontal alignment, it's pretty straight-forward. Italics
## might introduce some error (extending left of 0,0), but it's not a
## huge deal.

## For vertical alignment, there is a twist: the alignment can either
## be relative to the bounding box, or relative to the
## baseline. There are good use cases for both.

## Bounding box alignment is based on the text extents. If there are
## multiple lines, use the boundingRect(QRectF) method, otherwise
## tightBoundingRect. This is slow, but without it, the baseline
## effect would lead to misleading pictures.

## Baseline alignment uses only the ascent, and could be based on font
## metrics, or text extents (tightBoundingRect). R uses the font
## metrics, which is fastest. We will do the same.

## For multiple lines, each line is aligned the same as the block.

##' @param text The vector of strings to draw, recycled
##' @param halign The side of the text to horizontally align to the coordinate
##' @param valign The side of the text to vertically align to the
##' coordinate. Besides the obvious alignment options, there are two
##' different ways to center the text: according to the entire text
##' extents ("center") or only according to the region above the
##' baseline ("basecenter"). The former may be better for plotting
##' text, while the latter may be better for labeling.
##' @param rot The vector of the text rotations, in degrees, recycled
##' @param color The stroke color of the text
##' @param cex The vector of numeric expansion factors of the glyphs, recycled
##' @param hcex The vector of numeric horizontal expansion factors of
##' the glyphs, recycled. Overrides the \code{cex} in the horizontal
##' direction.
##' @param vcex The vector of numeric vertical expansion factors of
##' the glyphs, recycled. Overrides the \code{cex} in the vertical
##' direction.
##' @rdname painting 
qdrawText <- function(p, text, x, y, halign = c("center", "left", "right"),
                      valign = c("center", "basecenter", "baseline", "bottom",
                        "top"),
                      rot = 0, color = NULL, cex = 1.0, hcex = cex, vcex = cex)
{
  m <- max(length(text), length(x), length(y))
  text <- recycleVector(text, m)
  x <- recycleVector(x, m)
  y <- recycleVector(y, m)
  rot <- recycleVector(rot, m)
  hcex <- recycleVector(as.numeric(hcex), m)
  vcex <- recycleVector(as.numeric(vcex), m)
  drawText <- function(text, x, y, rot, color, hcex, vcex)
    invisible(.Call("qt_qdrawText_Painter", p, as.character(text),
                    as.numeric(x), as.numeric(y), as.integer(hflag + vflag),
                    as.numeric(rot), .normArgStroke(p, color, m), hcex, vcex))
  stopifnot(inherits(p, "Painter"))
  hflags <- c(left = 0x1, right = 0x2, center = 0x4)
  halign <- match.arg(halign)
  hflag <- hflags[halign]
  vflags <- c(top = 0x20, bottom = 0x40, center = 0x80)
  valign <- match.arg(valign)
  vflag <- vflags[valign]
  ## single lines should be vertically centered exactly
  if (valign == "center") {
    multi <- grepl("\n", text, fixed=TRUE)
    if (any(multi)) { ## draw the multilines immediately
      drawText(text[multi], x[multi], y[multi], rot[multi], color[multi],
               hcex[multi], vcex[multi])
      text <- text[!multi]; x <- x[!multi]; y <- y[!multi]
      rot <- rot[!multi]; color <- color[!multi];
      hcex <- hcex[!multi]; vcex <- vcex[!multi]
    }
    vflag <- NA
  }
  if (is.na(vflag)) {
    vflag <- vflags["top"]
    ascent <- qfontMetrics(p)["ascent"]
    adj <- 0
    if (valign == "basecenter")
      ascent <- ascent / 2
    else if (valign == "center") {
      extents <- qtextExtents(p, text)
      adj <- -(extents[,"y1"] - extents[,"y0"]) / 2
    }
    adj <- adj + ascent
    ## fix adjustment for rotation
    rads <- rot/360*2*pi
    tf <- qdeviceTransform(p)
    ## we perform an "inverse" rotation in Y, map Y to pixels, then back to X
    ## 'adj' is a magnitude, so we have to subtract the origin (0)
    ## this works around the flipped Y axis
    mapToX <- function(y)
      qvmap(tf$inverted(), qvmap(tf, 0, sin(rads) * y)[,2], 0)[,1]
    x <- x + mapToX(adj) - mapToX(0)
    y <- y + cos(rads)*adj
  }
  drawText(text, x, y, rot, color, hcex, vcex)
}

##' @rdname painting 
##' @param image A C++ \code{QImage} object
qdrawImage <- function(p, image, x, y) {
  stopifnot(inherits(p, "Painter"))
  stopifnot(inherits(image, "QImage"))
  m <- max(length(x), length(y))
  x <- recycleVector(x, m)
  y <- recycleVector(y, m)
  invisible(.Call("qt_qdrawImage_Painter", p, image, as.numeric(x),
                  as.numeric(y)))
}

##' @rdname painting 
qdrawGlyph <- function(p, path, x, y, cex = NULL, stroke = NULL, fill = NULL) {
  stopifnot(inherits(p, "Painter"))
  stopifnot(inherits(path, "QPainterPath"))
  m <- max(length(x), length(y))
  x <- recycleVector(x, m)
  y <- recycleVector(y, m)
  if (!is.null(cex)) {
    if (length(cex) == 1) {
      qglyphExpansion(p) <- cex
      cex <- NULL
    } else cex <- recycleVector(as.numeric(cex), m)
  }
  invisible(.Call("qt_qdrawGlyphs_Painter", p, path, as.numeric(x),
                  as.numeric(y), cex, .normArgStroke(p, stroke, m),
                  .normArgFill(p, fill, m)))
}

##' Get text extents and font metrics
##'
##' @title Text extents
##' @param p The paint context
##' @param text The text to analyze
##' @return A matrix representing the text bounds for
##' \code{qtextExtents}), a number for \code{qstrWidth} and
##' \code{qstrHeight}, or a vector with the ascent and descent for
##' \code{qfontMetrics}
##' @rdname text-extents
##' @author Michael Lawrence
qtextExtents <- function(p, text) {
  ans <- .Call("qt_qtextExtents_Painter", p, as.character(text))
  colnames(ans) <- c("x0", "y0", "x1", "y1")
  ans
}

##' @rdname text-extents
qstrWidth <- function(p, text) {
  ## FIXME: optimize by directly asking for widths, heights are expensive
  extents <- qtextExtents(p, text)
  extents[,3] - extents[,1]
}

##' @rdname text-extents
qstrHeight <- function(p, text) {
  extents <- qtextExtents(p, text)
  extents[,4] - extents[,2]
}

##' @rdname text-extents
qfontMetrics <- function(p) {
  stopifnot(inherits(p, "Painter"))
  ans <- .Call("qt_qfontMetrics_Painter", p)
  names(ans) <- c("ascent", "descent")
  ans
}
