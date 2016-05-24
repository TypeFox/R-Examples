#' Calculate colours from wavelength.
#'
#' \code{stat_color} computes color definitions according to human vision.
#'
#' @param mapping The aesthetic mapping, usually constructed with
#'    \code{\link[ggplot2]{aes}} or \code{\link[ggplot2]{aes_string}}. Only needs to be set
#'    at the layer level if you are overriding the plot defaults.
#' @param data A layer specific dataset - only needed if you want to override
#'    the plot defaults.
#' @param geom The geometric object to use display the data
#' @param position The position adjustment to use for overlapping points
#'    on this layer
#' @param show.legend logical. Should this layer be included in the legends?
#'   \code{NA}, the default, includes if any aesthetics are mapped.
#'   \code{FALSE} never includes, and \code{TRUE} always includes.
#' @param inherit.aes If \code{FALSE}, overrides the default aesthetics,
#'   rather than combining with them. This is most useful for helper functions
#'   that define both data and aesthetics and shouldn't inherit behaviour from
#'   the default plot specification, e.g. \code{\link[ggplot2]{borders}}.
#' @param ... other arguments passed on to \code{\link[ggplot2]{layer}}. This can
#'   include aesthetics whose values you want to set, not map. See
#'   \code{\link[ggplot2]{layer}} for more details.
#' @param na.rm	a logical value indicating whether NA values should be
#'   stripped before the computation proceeds.
#' @param type character one of "CMF" (color matching function) or "CC"
#'   (color coordinates).
#'
#' @section Computed variables:
#' \describe{
#'   \item{color}{color corresponding to x-value giving wavelength in
#'   nanometres.}
#' }
#'
#' @seealso \code{\link[photobiology]{color}}, which is used internally.
#'
#' @examples
#' library(photobiology)
#' library(ggplot2)
#' ggplot(sun.spct, aes(w.length, s.e.irrad)) + geom_line() +
#'   stat_color()
#'
#' @export
#' @family stats functions
#'
stat_color <- function(mapping = NULL, data = NULL, geom = "point",
                       type = "CMF",
                       position = "identity", na.rm = FALSE, show.legend = FALSE,
                       inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatColor, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(type = type,
                  na.rm = na.rm,
                  ...)
  )
}

#' @rdname gg2spectra-ggproto
#' @format NULL
#' @usage NULL
#' @export
#' @seealso \code{\link[ggplot2]{ggplot2-ggproto}}
StatColor <-
  ggplot2::ggproto("StatColor", ggplot2::Stat,
                   compute_group = function(data,
                                            scales,
                                            type) {
                   dplyr::mutate(data, color = photobiology::color(x, type))
                   },
                   default_aes = ggplot2::aes(color = ..color..,
                                              fill = ..color..),
                   required_aes = c("x", "y")
  )
