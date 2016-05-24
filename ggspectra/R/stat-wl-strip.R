#' Calculate colours from wavelength.
#'
#' \code{stat_wl_strip} computes color definitions according to human vision.
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
#' @param w.band photobiology::waveband object or a list of such objects or NULL.
#' @param length.out The number of steps to use to simulate a continuous
#'   range of colours when w.band == NULL.
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
#' # ggplot() methods for spectral objects set a default mapping for x and y.
#' ggplot(sun.spct) + geom_line() +
#'   stat_wl_strip(ymax = -0.02, ymin = -0.04) +
#'   scale_fill_identity()
#'
#' ggplot(sun.spct) + wl_guide(alpha = 0.33) + geom_line()
#'
#' @export
#' @family stats functions
#'
stat_wl_strip <- function(mapping = NULL, data = NULL, geom = "rect",
                       type = "CMF", w.band = NULL, length.out = 150,
                       position = "identity", na.rm = FALSE, show.legend = FALSE,
                       inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatColorGuide, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(type = type,
                  w.band = w.band,
                  length.out = length.out,
                  na.rm = na.rm,
                  ...)
  )
}

#' @rdname gg2spectra-ggproto
#' @format NULL
#' @usage NULL
#' @export
#' @seealso \code{\link[ggplot2]{ggplot2-ggproto}}
StatColorGuide <-
  ggplot2::ggproto("StatColorGuide", ggplot2::Stat,
                   compute_group = function(data,
                                            scales,
                                            type,
                                            w.band,
                                            length.out) {
                     if (is.null(w.band)) {
                       w.band <- split_bands(range(data$x), length.out = length.out)
                     } else {
                       w.band <- trim_waveband(w.band = w.band, range = data$x, trim = TRUE)
                     }

                     z <- wb2rect_spct(w.band = w.band)
                     names(z)[1] <- "x"
                     z
                    },
                   default_aes = ggplot2::aes(xmin = ..wl.low..,
                                              xmax = ..wl.high..,
                                              label = as.character(..wb.f..),
                                              fill = ..wb.color..),
                   required_aes = c("x")
  )

#' @rdname stat_wl_strip
#' @param ymin,ymax numeric used as aesthetics for plotting the guide.
#'
#' @export
#'
wl_guide <- function(mapping = NULL, data = NULL,
                             type = "CMF", w.band=NULL, length.out = 150,
                             ymin = -Inf, ymax = Inf,
                             position = "identity", na.rm = FALSE, show.legend = FALSE,
                             inherit.aes = TRUE, ...) {
  list(stat_wl_strip(mapping = mapping,
                        data = data,
                        geom = "rect",
                        type = type,
                        w.band = w.band,
                        length.out = length.out,
                        show.legend = show.legend,
                        inherit.aes = inherit.aes,
                        ymin = ymin,
                        ymax = ymax,
                        ...),
       scale_fill_identity()
  )
}
