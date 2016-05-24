## convenience wrapper for constructor

##' Create a scene, a container of layers
##'
##' Often, a \code{PlotView} object is passed as \code{parent}, so
##' that the scene is not deleted until the view is destroyed.
##' @title Create scenes
##' @param parent If non-NULL, a C++ \code{QObject} object that owns
##' the memory of the scene
##' @return A C++ \code{QGraphicsScene} object
##' @author Michael Lawrence
qscene <- function(parent = NULL) Qt$QGraphicsScene(parent)

##' These functions constitute a vectorized API for adding graphical
##' primitive items into a scene. Each primitive is represented by a
##' separate canvas item. This is a different paradigm from the
##' layer-based painting, which has a more appropriate design for most
##' interactive graphics. Most of the time, the user should ignore this API.
##'
##' 
##' @title Adding scene items
##' @param s The scene, a C++ \code{QGraphicsScene}
##' @param x The X coordinates, as expected by \code{\link{xy.coords}}, recycled
##' @param y The Y coordinates, as expected by \code{\link{xy.coords}}, recycled
##' @param radius The scalar radius for the points, in device coordinates
##' @author Deepayan Sarkar
##' @rdname add-scene-items
qscene.points <- function(s, x, y, radius = 1)
{
  xy <- xy.coords(x, y, recycle = TRUE)
  .Call("scene_addPoints", s, as.double(xy$x), as.double(xy$y),
        as.double(radius), PACKAGE = "qtpaint")
}

##' @rdname add-scene-items
##' @param lwd The vector of line widths, in device coordinates, recycled
qscene.lines <- function(s, x, y, lwd = 0)
{
  xy <- xy.coords(x, y, recycle = TRUE)
  .Call("scene_addLines", s, as.double(xy$x), as.double(xy$y), as.double(lwd),
        PACKAGE = "qtpaint")
}

##' @param x1 The first X coordinates, recycled
##' @param y1 The first Y coordinates, recycled
##' @param x2 The second X coordinates, recycled
##' @param y2 The second Y coordinates, recycled
##' @rdname add-scene-items
qscene.segments <- function(s, x1, y1, x2, y2, lwd = 0)
{
  n <- max(length(x1), length(x2), length(y1), length(y2))
  .Call("scene_addSegments",
        s,
        rep(as.double(x1), length.out = n),
        rep(as.double(y1), length.out = n),
        rep(as.double(x2), length.out = n),
        rep(as.double(y2), length.out = n),
        as.double(lwd),
        PACKAGE = "qtpaint")
}

##' @param w Vector of rectangle widths, recycled
##' @param h Vector of rectangle heights, recycled
##' @rdname add-scene-items
qscene.rect <- function(s, x, y, w = 1, h = 1)
{
  xy <- xy.coords(x, y, recycle = TRUE)
  .Call("scene_addRect", s, as.double(xy$x), as.double(xy$y), as.double(w),
        as.double(h), PACKAGE = "qtpaint")
}

##' Set properties of all items in a scene. Fast path for when a scene
##' contains many items.
##'
##' @title Item properties
##' @param x The scene, a C++ \code{QGraphicsScene}
##' @param flag A value or combination of values from the
##' \code{QGraphicsItem::GraphicsItemFlag} enumeration
##' @param status Whether the flag should be set to \code{TRUE} or \code{FALSE}
##' @author Deepayan Sarkar and Michael Lawrence
##' @rdname item-properties
qsetItemFlags <- function(x, flag = Qt$QGraphicsItem$ItemIsMovable,
                          status = FALSE)
{
  .Call("qt_qsetItemFlags", x, flag, status, PACKAGE = "qtpaint")
}

##' @param mode Whether a text item (C++ \code{QGraphicsTextItem})
##' should behave as a text editor, text browser, or not allow
##' interaction.
##' @rdname item-properties
qsetTextItemInteraction <- function(x, mode = c("none", "editor", "browser"))
{
  mode <- match.arg(mode)
  .Call("qt_qsetTextItemInteraction", x, mode, PACKAGE = "qtpaint")
}
