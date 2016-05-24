## Conveniences for interacting with Qanviz::Layer.

##' This function constructs a Layer, an item in the canvas that may
##' paint part of the plot and/or respond to user input. The behavior
##' of the layer is implemented using R functions, which are passed to
##' the constructor. Other arguments determine the scaling and
##' positioning of the layer, and whether the drawing of the layer is
##' cached and/or clipped.
##'
##' All drawing and user input handling is performed by R callbacks,
##' which must accept a specific set of arguments. The paint callback,
##' passed as \code{paintFun}, must take at least two arguments,
##' conventionally named \code{layer} and \code{painter}. The
##' \code{layer} argument is a C++ \code{RLayer} object, the same
##' instance that was created by calling the constructor. All painting
##' is performed through the \code{painter} argument, which is a C++
##' \code{Painter} object. See the
##' \code{\link[qdrawGlyph]{paint functions}} for more details. The
##' \code{paintFun} may take one additional, optional argument,
##' conventionally named \code{exposed}, which is the rectangle, in
##' layer coordinates, that needs to be drawn.
##'
##' All of the other callbacks, except for \code{sizeHintFun}, are
##' event handlers. Two arguments are passed, conventionally named
##' \code{layer} and \code{event}. The former is the layer constructed
##' in the call to \code{qlayer}, and the latter describes the event
##' as an instance of C++ \code{QGraphicsSceneEvent}. The exact
##' subclass depends on the event. Manipulating an event currently
##' requires low-level calls through the \code{qtbase} package. See
##' its documentation.
##' 
##' @title Create a layer
##' @param parent The scene, for a top-level layer, or the parent
##' layer that contains the new layer in a grid layout
##' @param paintFun The function that implements painting, called
##' whenever the layer needs to be repainted. All drawing occurs in
##' data/layer coordinates.
##' @param keyPressFun The function called when a key is pressed
##' @param keyReleaseFun The function called when a key is released
##' @param mouseDoubleClickFun The function called when a mouse button
##' is double-clicked
##' @param mouseMoveFun The function called when the mouse is moved
##' while holding down a button 
##' @param mousePressFun The function called when a mouse button is pressed
##' @param mouseReleaseFun The function called when a mouse button is released
##' @param wheelFun The function called when the mouse wheel is turned
##' @param hoverMoveFun The function called when the mouse moves
##' without any buttons pressed
##' @param hoverEnterFun The function called when the mouse pointer
##' enters the layer
##' @param hoverLeaveFun The function called when the mouse pointer
##' leaves the layer
##' @param contextMenuFun The function called when a context menu is
##' requested, such as through a right mouse button click
##' @param dragEnterFun The function called when the pointer enters
##' the layer while dragging something
##' @param dragLeaveFun The function called when the pointer leaves
##' the layer while dragging something
##' @param dragMoveFun The function called when the pointer moves within
##' the layer while dragging something
##' @param dropFun The function called when something is dropped on the layer
##' @param focusInFun The function called when the layer gains the
##' keyboard focus
##' @param focusOutFun The function called when the layer loses the
##' keyboard focus
##' @param sizeHintFun The function called to determine the size
##' constraints of the layer. Not currently documented.
##' @param limits A \code{QRectF}, possibly created by
##' \code{\link[qtbase]{qrect}}, indicating the X and Y scales of the
##' layer in data/layer coordinates
##' @param row The 0-based row index of the layer in the parent grid layout
##' @param col The 0-based column index of the layer in the parent grid layout
##' @param rowSpan The 0-based number of rows spanned by the layer in the layout
##' @param colSpan The 0-based number of cols spanned by the layer in the layout
##' @param geometry A \code{QRectF}, possibly created by
##' \code{\link[qtbase]{qrect}}, indicating the position and size of
##' the layer in figure/scene coordinates. This is overridden by the
##' parent grid layout, so is really only useful for a top-level
##' layer. A warning will be issued if the geometry is specified along
##' with a parent layer. We also issue a warning if this argument is
##' specified when the scene has a view in "geometry" rescale mode,
##' because the view determines the geometry. The default geometry is
##' the bounding rectangle of the scene, if not null, or 600x400
##' otherwise.
##' @param clip Logical indicating whether to clip drawing to within the layer
##' @param cache Logical indicating whether to cache drawing, which
##' helps performance for relatively static layers sitting under more
##' dynamic ones
##' @return The layer, a C++ instance of \code{RLayer}
##' @author Michael Lawrence
##' @examples
##' scene <- qscene()
##' layer <- qlayer(scene, function(layer, painter) {
##'   qdrawCircle(painter, 1:10, 1:10, 1)
##' }, limits = qtbase::qrect(0, 0, 11, 11))
##' qplotView(scene)
qlayer <- function(parent = NULL, paintFun = NULL, keyPressFun = NULL,
                   keyReleaseFun = NULL, mouseDoubleClickFun = NULL,
                   mouseMoveFun = NULL, mousePressFun = NULL,
                   mouseReleaseFun = NULL, wheelFun = NULL,
                   hoverMoveFun = NULL, hoverEnterFun = NULL,
                   hoverLeaveFun = NULL, contextMenuFun = NULL,
                   dragEnterFun = NULL, dragLeaveFun = NULL,
                   dragMoveFun = NULL, dropFun = NULL,
                   focusInFun = NULL, focusOutFun = NULL,
                   sizeHintFun = NULL, limits = qrect(),
                   row = 0L, col = 0L, rowSpan = 1L, colSpan = 1L,
                   geometry = defaultLayerGeometry(parent), clip = cache,
                   cache = FALSE)
{
  if (cache && !clip)
    warning("Enabling caching implicitly enables clipping")
  p <- NULL
  if (inherits(parent, "QGraphicsItem"))
    p <- parent
  args <- list(p,
               .normArgCallback(paintFun),
               .normArgCallback(keyPressFun),
               .normArgCallback(keyReleaseFun),
               .normArgCallback(mouseDoubleClickFun),
               .normArgCallback(mouseMoveFun),
               .normArgCallback(mousePressFun),
               .normArgCallback(mouseReleaseFun),
               .normArgCallback(wheelFun),
               .normArgCallback(hoverMoveFun),
               .normArgCallback(hoverEnterFun),
               .normArgCallback(hoverLeaveFun),
               .normArgCallback(contextMenuFun),
               .normArgCallback(dragEnterFun),
               .normArgCallback(dragLeaveFun),
               .normArgCallback(dragMoveFun),
               .normArgCallback(dropFun),
               .normArgCallback(focusInFun),
               .normArgCallback(focusOutFun),               
               .normArgCallback(sizeHintFun))
  layer <- .Call("qanviz_RLayer", args, PACKAGE="qtpaint")
  if (inherits(parent, "Qanviz::Layer")) {
    if (!missing(geometry)) {
      warning("geometry will be overridden by parent layout")
    }
    parent$addLayer(layer, row, col, rowSpan, colSpan)
  } else if (inherits(parent, "QGraphicsScene")) {
    parent$addItem(layer)
    viewGeometry <- viewGeometry(parent)
    if (!is.null(viewGeometry)) {
      if (!missing(geometry)) {
        warning("geometry will be overridden by view in geometry rescale mode")
      }
      geometry <- viewGeometry
    }
  } else if (!is.null(parent)) stop("Unsupported parent type")
  layer$geometry <- geometry
  layer$setLimits(limits)
  layer$setFlag(Qt$QGraphicsItem$ItemClipsToShape, clip)
  if (!cache)
    layer$setCacheMode(Qt$QGraphicsItem$NoCache)
  layer
}

defaultLayerGeometry <- function(parent) {
  geometry <- qrect(0, 0, 600, 400)
  if (is(parent, "QGraphicsScene")) {
    if (!parent$sceneRect$isNull()) {
      geometry <- parent$sceneRect
    }
  }
  geometry
}

viewGeometry <- function(scene) {
  views <- scene$views()
  views <- Filter(function(v) v$rescaleMode() == Qanviz$PlotView$WidgetGeometry,
                  views)
  if (length(views) > 0L) {
    views[[1L]]$viewport()$rect
  }
}

.normArgCallback <- function(callback) {
  if (!is.null(callback))
    callback <- as.function(callback)
  callback
}

##' Add or retrieve a layer to or from the grid layout of another layer 
##'
##' @title Grid layout accessors
##' @param x Parent layer
##' @param i 0-based row position
##' @param j 0-based column position
##' @param rowSpan Number of rows spanned by \code{value}
##' @param colSpan Number of columns spanned by \code{value}
##' @param value The layer to add to the layout
##' @author Michael Lawrence
##' @rdname layout-accessors
##' @method [<- `Qanviz::Layer`
"[<-.Qanviz::Layer" <- function (x, i = 0, j = 0, rowSpan = 1, colSpan = 1,
                                 value)
{
  x$addLayer(value, i, j, rowSpan, colSpan)
  x
}
##' @rdname layout-accessors
##' @method [ `Qanviz::Layer`
"[.Qanviz::Layer" <-
  function (x, i = 0, j = 0)
{
  x$layerAt(i, j)
}

## wrappers that could be added by anyone interested in maintaining them:

## qlimits <- function(x) x$limits()
## "qlimits<-" <- function(x, limits) x$setLimits(limits)
## qaddLayer <- function(x, child) x$addLayer(child)
## qdeviceTransform: layer$deviceTransform(event)
## q(row,col)Stretch[<-]: layer$layout()$(row,col)[Set]Stretch()
## q(h,v)Spacing[<-]: layer$layout()$setHorizontal(Vertical)Spacing()
## qbackgroundBrush[<-]: scene$backgroundBrush <- brush
## qclearSelection: layer$clearSelection()
## qzValue[<-]: layer$setZValue()
## qlocate: layer$locate()
## qminimumSize[<-]: layer$setMinimumSize()
## qcacheMode[<-]: layer$setCacheMode()
## qclip[<-]: layer$setClip()
## qfocus[<-]: layer$setFocus()
## qoverlayScene: layer$overlayScene()
