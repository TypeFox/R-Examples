#' Integrate ranges under spectral curve.
#'
#' \code{stat_wb_contribution} computes means under a curve. It first integrates the
#'   area under a spectral curve and also the mean expressed per nanaometre of
#'   wavelength for each waveband in the input. Sets suitable default aestheics
#'   for "rect", "hline", "vline", "text" and "label" geoms displaying "contributions"
#'   per waveband to the total of the spectral integral.
#'
#' @param mapping The aesthetic mapping, usually constructed with
#'   \code{\link[ggplot2]{aes}} or \code{\link[ggplot2]{aes_string}}. Only needs
#'   to be set at the layer level if you are overriding the plot defaults.
#' @param data A layer specific dataset - only needed if you want to override
#'   the plot defaults.
#' @param geom The geometric object to use display the data
#' @param position The position adjustment to use for overlapping points on this
#'   layer
#' @param show.legend logical. Should this layer be included in the legends?
#'   \code{NA}, the default, includes if any aesthetics are mapped. \code{FALSE}
#'   never includes, and \code{TRUE} always includes.
#' @param inherit.aes If \code{FALSE}, overrides the default aesthetics, rather
#'   than combining with them. This is most useful for helper functions that
#'   define both data and aesthetics and shouldn't inherit behaviour from the
#'   default plot specification, e.g. \code{\link[ggplot2]{borders}}.
#' @param ... other arguments passed on to \code{\link[ggplot2]{layer}}. This
#'   can include aesthetics whose values you want to set, not map. See
#'   \code{\link[ggplot2]{layer}} for more details.
#' @param na.rm	a logical value indicating whether NA values should be stripped
#'   before the computation proceeds.
#' @param w.band a waveband object or a list of waveband objects or numeric
#'   vector of at least length two.
#' @param integral.fun function on $x$ and $y$.
#' @param label.mult numeric Scaling factor applied to y-integral values before
#'   conversion into character strings.
#' @param label.fmt character string giving a format definition for converting
#'   y-integral values into character strings by means of function
#'   \code{\link{sprintf}}.
#' @param ypos.mult numeric Multiplier constant used to scale returned
#'   \code{y} values.
#' @param ypos.fixed numeric If not \code{NULL} used a constant value returned
#'   in \code{y}.
#'
#' @section Computed variables:
#' \describe{
#'   \item{label}{intergral value as formatted text}
#'   \item{x}{w.band-midpoint}
#'   \item{xmin}{w.band minimum}
#'   \item{xmax}{w.band maximum}
#'   \item{ymean}{Mean value as numeric}
#'   \item{yint}{Integral value as numeric}
#' }
#'
#' @import photobiology
#'
#' @examples
#' library(photobiology)
#' library(photobiologyWavebands)
#' library(ggplot2)
#' # ggplot() methods for spectral objects set a default mapping for x and y.
#' ggplot(sun.spct) + geom_line() +
#'   stat_wb_contribution(w.band = VIS_bands())
#'
#' @export
#' @family stats functions
#'
stat_wb_contribution <- function(mapping = NULL, data = NULL, geom = "rect",
                       w.band = NULL,
                       integral.fun = photobiology::integrate_xy,
                       label.mult = 1,
                       label.fmt = "%.3g",
                       ypos.mult = 1.07,
                       ypos.fixed = NULL,
                       position = "identity", na.rm = FALSE, show.legend = NA,
                       inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatWbContrib, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(w.band = w.band,
                  integral.fun = integral.fun,
                  label.mult = label.mult,
                  label.fmt = label.fmt,
                  ypos.mult = ypos.mult,
                  ypos.fixed = ypos.fixed,
                  na.rm = na.rm,
                  ...)
  )
}

#' @rdname gg2spectra-ggproto
#' @format NULL
#' @usage NULL
#' @export
StatWbContrib <-
  ggplot2::ggproto("StatWbContrib", ggplot2::Stat,
                   compute_group = function(data,
                                            scales,
                                            w.band,
                                            integral.fun,
                                            label.mult,
                                            label.fmt,
                                            ypos.mult,
                                            ypos.fixed) {
                     if (is.null(w.band)) {
                       w.band <- waveband(data$x)
                     }
                     if (is.any_spct(w.band) ||
                         (is.numeric(w.band) && length(na.omit(w.band)) >= 2)) {
                       w.band <- waveband(range(w.band, na.rm = TRUE))
                     }
                     if (!is.list(w.band) || is.waveband(w.band)) {
                       w.band <- list(w.band)
                     }
                     stopifnot(is.function(integral.fun))
                     w.band <- trim_wl(w.band, data$x)
                     integ.df <- data.frame()
                     for (wb in w.band) {
                       if (is.numeric(wb)) { # user supplied a list of numeric vectors
                         wb <- waveband(wb)
                       }

                       range <- range(wb)
                       mydata <- trim_tails(data$x, data$y, use.hinges = TRUE,
                                            low.limit = range[1],
                                            high.limit = range[2])
                       if (is_effective(wb)) {
                         warning("BSWFs not supported by summary: using wavelength range for ",
                                 labels(wb)$label, "'.")
                         wb <- waveband(wb)
                       }
                       yint.tmp <- integral.fun(mydata$x, mydata$y)
                       ymean.tmp <- yint.tmp / spread(wb)
                       integ.df <- rbind(integ.df,
                                         data.frame(x = midpoint(mydata$x),
                                                    xmin = min(wb),
                                                    xmax = max(wb),
                                                    ymin = min(data$y),
                                                    ymax = max(data$y),
                                                    yint = yint.tmp,
                                                    ymean = ymean.tmp,
                                                    wb.color = color(wb),
                                                    wb.name = labels(wb)$label)
                                         )
                     }
                     if (is.null(ypos.fixed)) {
                       integ.df$y <- with(integ.df, ymin + (ymax - ymin) * ypos.mult)
                     } else {
                       integ.df$y <- ypos.fixed
                     }
                     integ.df$yint <- integ.df$yint / integral.fun(data$x, data$y)
                     integ.df$y.label <- sprintf(label.fmt, integ.df$yint * label.mult)
#                     print(integ.df)
                     integ.df
                   },
                   default_aes = ggplot2::aes(label = ..y.label..,
                                              xmin = ..xmin..,
                                              xmax = ..xmax..,
                                              ymin = ..y.. - (..ymax.. - ..ymin..) * 0.03,
                                              ymax = ..y.. + (..ymax.. - ..ymin..) * 0.03,
                                              yintercept = ..ymean..,
                                              fill = ..wb.color..),
                   required_aes = c("x", "y")
  )
