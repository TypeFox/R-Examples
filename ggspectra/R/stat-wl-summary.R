#' Average area under curve for regions.
#'
#' \code{stat_wl_summary} computes the area under a curve.
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
#' @param range a numeric vector of at least length two.
#' @param integral.fun function on $x$ and $y$.
#' @param label.fmt character string giving a format definition for converting
#'   y-integral values into character strings by means of function \code{\link{sprintf}}.
#' @section Computed variables:
#' \describe{
#'   \item{label}{intergral value as formatted text}
#'   \item{x}{range-midpoint}
#'   \item{xmin}{range minimum}
#'   \item{xmax}{range maximum}
#'   \item{y}{integral value as numeric}
#' }
#'
#' @examples
#' library(photobiology)
#' library(ggplot2)
#' # ggplot() methods for spectral objects set a default mapping for x and y.
#' ggplot(sun.spct) + geom_line() + stat_wl_summary(geom = "hline")
#' ggplot(sun.spct) + geom_line() +
#'  stat_wl_summary(label.fmt = "%.3f", color = "red")
#'
#' @export
#' @family stats functions
#'
stat_wl_summary <- function(mapping = NULL, data = NULL, geom = "text",
                       range = NULL,
                       integral.fun = photobiology::integrate_xy, label.fmt = "%.3g",
                       position = "identity", na.rm = FALSE, show.legend = NA,
                       inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatAverage, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(range = range,
                  integral.fun = integral.fun,
                  label.fmt = label.fmt,
                  na.rm = na.rm,
                  ...)
  )
}

#' @rdname gg2spectra-ggproto
#' @format NULL
#' @usage NULL
#' @export
StatAverage <-
  ggplot2::ggproto("StatAverage", ggplot2::Stat,
                   compute_group = function(data,
                                            scales,
                                            range,
                                            integral.fun,
                                            label.fmt,
                                            summary.fmt) {
                     range <-
                       photobiology::normalize_range_arg(range,
                                                         default.range = range(data$x))
                     mydata <- photobiology::trim_tails(data$x, data$y,
                                                        low.limit = range[1],
                                                        high.limit = range[2])
                     integ.df <- data.frame(x = midpoint(mydata$x),
                                            xmin = range[1],
                                            xmax = range[2],
                                            y = integral.fun(mydata$x, mydata$y) /
                                              (range[2] - range[1]))
                     integ.df$label <- sprintf(label.fmt, integ.df$y)
                     integ.df
                   },
                   default_aes = ggplot2::aes(label = ..label..,
                                              x = ..x..,
                                              xmin = ..xmin..,
                                              xmax = ..xmax..,
                                              y = ..y..,
                                              ymax = ..y..,
                                              ymin = 0 * ..y..,
                                              yintercept = ..y..),
                   required_aes = c("x", "y")
  )

