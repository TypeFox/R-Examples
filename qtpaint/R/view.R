##' An extension of \code{QGraphicsView} with special functionality
##' for plotting statistical graphics. 
##'
##' The \code{PlotView} class extends \code{QGraphicsView} to add two
##' new features, from the user perspective. First, it automates
##' rescaling when the widget is resized. There are three rescale
##' modes: \code{geometry}, \code{transform} and \code{none}:
##' 
##' \describe{
##' \item{\code{geometry}}{Most similar to how base R
##' graphics devices behave: the geometry of the figure is fixed to
##' the geometry of the view. This mode is often convenient but really
##' only works if there is only a single view for each scene.}
##' 
##' \item{\code{transform}}{Adjusts the view transform so that the
##' visible region of the scene does not change due to size
##' changes. This is most similar to the behavior of GGobi and
##' supports multiple views of the same scene. The downside is that
##' the layout is not activated, so it cannot adapt to make better use
##' of the available space. It works best with graphics that overlay
##' guides on the plot, rather than position them adjacent to the plot
##' in a layout.}
##'
##' \item{\code{none}}{No rescaling is performed; when the size
##' changes, more or less of the plot is shown. This is probably the most
##' common mode in zoomable user interfaces.}
##' }
##'
##' The other feature is the overlay scene: a separate scene that is
##' fixed to the geometry of the viewport. It is always shown over the
##' primary scene and it is stationary across transformations
##' and scrolling of the viewport. This is useful for overlaying
##' guides on a plot in a fixed position, like the axes in GGobi. Call
##' the \code{overlay} method on a plot view instance to obtain the
##' overlay scene and manipulate it directly.
##' 
##' @title Plot view
##' @param scene The scene, a \code{QGraphicsScene}
##' @param parent The parent \code{QObject}, usually a
##' \code{QWidget} to contain the view
##' @param rescale The rescale mode, see details
##' @param opengl If \code{TRUE}, use OpenGL, otherwise the software driver
##' @return A C++ \code{PlotView} object
##' @author Michael Lawrence
qplotView <- function(scene, parent = NULL,
                      rescale = c("geometry", "transform", "none"),
                      opengl = TRUE)
{
  rescale <- c(none = 0L, geometry = 1L, transform = 2L)[match.arg(rescale)]
  view <- Qanviz$PlotView(scene, parent, rescale, opengl)
  if (is.null(scene$parent())) # view becomes default parent of scene
    scene$setParent(view)
  view
}
