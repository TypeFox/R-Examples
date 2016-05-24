# GeomAxis <- proto(ggplot2:::Geom, {
#   draw <- function(., data, scales, coordinates, location = c(0.2, 0.2), size = 0.9, colour = "black", axis, ...) {
#     axesVp <- axesViewport(axis, location)
#     axes <- axesGrob(axis, gp = gpar(col = colour))
#
#     gTree(
#       children = gList(axes),
#       childrenvp = axesVp
#     )
#   }
#
#   objname <- "axis"
#   desc <- "Projection axes"
#
#   default_stat <- function(.) StatIdentity
#   required_aes <- c()
#   default_aes <- function(.) aes()
#
# })

#' Geom Axis.
#'
#' A special ggplot2 geom for drawing the tour axes
#'
#' @param ... should include data, location, aes\_string information
#' @author Hadley Wickham \email{h.wickham@@gmail.com}
#' @keywords internal
#' @export
#' @examples
#' library(ggplot2)
#' print(ggplot(dd_example("tour2d")))
geom_axis <- function(axis, location, ...) {

  GeomTourAxis <- ggproto("GeomTourAxis", Geom,
    required_aes = c(),
    default_aes = aes(),
    draw_key = NULL,

    draw_panel = function(ignore_data, ignore_panel_scales, ignore_coord) {

      # axis and location will be filled b/c of scope,
      # even though they're not supplied
      axesVp <- axesViewport(axis, location)
      axes <- axesGrob(axis, gp = gpar(col = axis[, "colour"]))

      gTree(
        children = gList(axes),
        childrenvp = axesVp
      )
    }
  )


  layer(
    geom = GeomTourAxis, data = axis, inherit.aes = FALSE,
    # remaining args to make the layer happy
    params = list(na.rm = FALSE), mapping = aes(),  stat = "identity",
    position = "identity", show.legend = FALSE
  )
}
