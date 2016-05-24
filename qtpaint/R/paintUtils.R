### some utilities to help make drawing easier

##' Some glyph constructors for simple glyphs like circle, square and
##' triangle. For use with \code{\link{qdrawGlyph}}.
##'
##' @title Glyph constructors
##' @param r Radius of the circle
##' @return C++ \code{QPainterPath} instance for use with
##' \code{\link{qdrawGlyph}}
##' @author Michael Lawrence
##' @rdname qglyphs
qglyphCircle <- function(r = 5) {
  glyph <- Qt$QPainterPath()
  glyph$addEllipse(qpoint(0, 0), r, r)
  glyph
}
##' @param x Length of one side of the square or triangle or the segment
##' @rdname qglyphs
qglyphSquare <- function(x = 5) {
  glyph <- Qt$QPainterPath()
  glyph$addRect(-x / 2, -x / 2, x, x)
  glyph
}
##' @param direction Whether the triangle is pointing up or down, or
##' the direction of the segment (unit: radian)
##' @rdname qglyphs
qglyphTriangle <- function(x = 5, direction = c("up", "down")) {
  direction <- match.arg(direction)
  if (direction == "down")
    x <- -x
  glyph <- Qt$QPainterPath()
  glyph$moveTo(-x, x)
  glyph$lineTo(x, x)
  glyph$lineTo(0, -x)
  glyph$closeSubpath()
  glyph
}
##' @param text The text of the text glyph
##' @param size The font size of the text glyph
##' @rdname qglyphs
qglyphText <- function(text = "X", size = 12) {
  glyph <- Qt$QPainterPath()
  glyph$addText(-size / 2, size / 2, qfont(pointsize = size), text)
  glyph
}
##' @rdname qglyphs
qglyphSegment <- function(x = 5, direction = 0) {
  glyph <- Qt$QPainterPath()
  x0 <- x * cos(direction)
  y0 <- x * sin(direction)
  glyph$moveTo(-x0, -y0)
  glyph$lineTo(x0, y0)
  glyph
}

##' Transforms X and Y coordinates with a Qt-style transformation
##' matrix. The advantage over direct use of Qt is vectorization.
##'
##' @title Mapping coordinates
##' @param m A matrix encoding the transformation, or something
##' coercible to a matrix, like a C++ \code{QTransform} instance
##' @param x X coordinates; if \code{y} is missing, should be
##' something coercible to a numeric vector or matrix. If the vector
##' coercion succeeds, the vector is coerced to a matrix with
##' \code{matrix(x, ncol = 2, byrow = TRUE)}. The first column is taken
##' as X, the second as Y.
##' @param y Y coordinates, optional
##' @return The mapped coordinates, as a two column (X, Y) matrix,
##' unless \code{y} is missing, in which case an attempt is made to
##' coerce the result to the class of \code{x}, if any.
##' @author Michael Lawrence
qvmap <- function(m, x, y) {
  m <- as.matrix(m)
  cl <- NULL
  if (missing(y)) {
    cl <- class(x) # if only 'x' specified, preserve its class
    if (NCOL(x) == 1L) {
      tmpX <- try(as.numeric(x), silent=TRUE)
      if (!is(tmpX, "try-error"))
        x <- matrix(tmpX, ncol = 2, byrow = TRUE)
      else {
        x <- try(as.matrix(x), silent=TRUE)
        if (is(x, "try-error"))
          stop("if 'y' is missing, 'x' must be coercible to vector or matrix")
      }
    }
    y <- x[,2]
    x <- x[,1]
  }
  mapped <- cbind(x * m[1,1] + y * m[2,1] + m[3,1],
                  y * m[2,2] + x * m[1,2] + m[3,2])
  if (!is.null(cl) && canCoerce(mapped, cl))
    mapped <- as(mapped, cl)
  mapped
}

## Creates a QTransform that flips the Y axis

.validRect <- function(r) {
  is.matrix(r) && is.numeric(r) && identical(dim(r), c(2L, 2L))
}

##' Generate transform for flipping Y axis.
##'
##' @title Flip the Y axis
##' @param ymax Maximum Y value or a rectangle (\code{QRect} or matrix)
##' @param ymin Minimum Y value
##' @return A \code{QTransform} object that will transform points by
##' flipping the axis.
##' @seealso \code{\link{qvmap}}
##' @author Michael Lawrence
##' @rdname qflipy
qflipY <- function(ymax, ymin = 0) UseMethod("qflipY")
##' @method qflipY numeric
##' @rdname qflipy
qflipY.numeric <- function(ymax, ymin = 0) {
  if (.validRect(ymax)) {
    ymin <- ymax[3]
    ymax <- ymax[4]
  }
  Qt$QTransform(1, 0, 0, -1, 0, (ymax + ymin))
}
##' @method qflipY QRect
##' @rdname qflipy
qflipY.QRectF <- function(ymax, ymin = 0) qflipY(as.matrix(ymax))
##' @method qflipY QRectF
##' @rdname qflipy
qflipY.QRect <- qflipY.QRectF

##' Force a redraw of a layer or scene, clearing the cache. This needs
##' to be called whenever the drawing would change, e.g., if the data
##' or some visual attribute has changed. There is no automatic way
##' for qtpaint to detect this.
##'
##' @title Updating drawings
##' @param x The object, usually a layer or scene, to be redrawn
##' @author Michael Lawrence
##' @rdname qupdate
qupdate <- function(x) UseMethod("qupdate")

##' @method qupdate QGraphicsView
##' @rdname qupdate
qupdate.QGraphicsView <- function(x) {
  qupdate(x$scene())
  x$viewport()$repaint()
}

refreshItemCache <- function(item) {
  mode <- item$cacheMode()
  item$setCacheMode(Qt$QGraphicsItem$NoCache)
  item$setCacheMode(mode)
}

##' @method qupdate QGraphicsScene
##' @rdname qupdate
qupdate.QGraphicsScene <- function(x) {
  lapply(x$items(), refreshItemCache)
  x$update()
}

##' @method qupdate QGraphicsItem
##' @rdname qupdate
qupdate.QGraphicsItem <- function(x) {
  refreshItemCache(x)
  x$update()
}
