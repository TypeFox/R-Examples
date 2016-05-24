#' Create a new ggplot plot from spectral data.
#'
#' \code{ggplot()} initializes a ggplot object. It can be used to
#' declare the input spectral object for a graphic and to optionally specify the
#' set of plot aesthetics intended to be common throughout all
#' subsequent layers unless specifically overridden.
#'
#' \code{ggplot()} is typically used to construct a plot
#' incrementally, using the + operator to add layers to the
#' existing ggplot object. This is advantageous in that the
#' code is explicit about which layers are added and the order
#' in which they are added. For complex graphics with multiple
#' layers, initialization with \code{ggplot} is recommended.
#'
#' There are three common ways to invoke \code{ggplot}:
#' \itemize{
#'    \item \code{ggplot(spct, aes(x, y, <other aesthetics>))}
#'    \item \code{ggplot(spct, unit.out = <unit.to.use>)}
#'    \item \code{ggplot(spct)}
#'   }
#' The first method is recommended if all layers use the same
#' data and the same set of aesthetics, although this method
#' can also be used to add a layer using data from another
#' data frame. See the first example below. The second
#' method specifies the default spectrum object to use for the plot, and the
#' units to be used for y in the plot,
#' but no aesthetics are defined up front. This is useful when
#' one data frame is used predominantly as layers are added,
#' but the aesthetics may vary from one layer to another. The
#' third method specifies the default spectrum object to use for the plot,
#' but no aesthetics are defined up front. This is useful when
#' one spectrum is used predominantly as layers are added,
#' but the aesthetics may vary from one layer to another.
#'
#' @param data Default spectrum dataset to use for plot. If not a spectrum, the
#'   methods used will be those defined in package \code{ggplot2}. See \code{\link[ggplot2]{ggplot}}.
#'   If not specified,
#'   must be suppled in each layer added to the plot.
#' @param mapping Default list of aesthetic mappings to use for plot.
#'   If not specified, in the case of spectral objects, a default mapping will
#'   be used.
#' @param unit.out character string indicating type of units to use for
#'   plotting.
#' @param ... Other arguments passed on to methods. Not currently used.
#' @param environment If an variable defined in the aesthetic mapping is not
#'   found in the data, ggplot will look for it in this environment. It defaults
#'   to using the environment in which \code{ggplot()} is called.
#' @export
#' @examples
#' library(photobiology)
#' library(ggplot2)
#' ggplot(sun.spct) + geom_line()
#' ggplot(sun.spct, unit.out = "photon") + geom_line()
#'
#' ggplot(yellow_gel.spct) + geom_line()
#' ggplot(yellow_gel.spct, plot.qty = "absorbance") + geom_line()
#'
#' @note Current implementation does not merge default mapping with user
#' supplied mapping. If user supplies a mapping, it is used as is, and
#' variables should be present in the spectral object. In contrast, when
#' using the default mapping, unit conversion is done on the fly when needed.
#' To add to the default mapping, aes() can be used by itself to compose
#' the ggplot.
#'
#' @name ggplot
#'
ggplot.source_spct <-
  function(data, mapping = NULL, ...,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           environment = parent.frame()) {
    if (is.null(mapping)) {
    if (unit.out == "energy") {
      data <- q2e(data)
      mapping <- aes_(~w.length, ~s.e.irrad)
    } else if (unit.out %in% c("photon", "quantum")) {
      data <- e2q(data)
      mapping <- aes_(~w.length, ~s.q.irrad)
    } else {
      stop("Invalid 'unit.out' argument value: '", unit.out, "'")
    }
  }
  rmDerivedSpct(data)
  ggplot(data = data, mapping = mapping, ...,
                    environment = environment)
}
#' @rdname ggplot
#'
#' @export
#'
ggplot.response_spct <-
  function(data, mapping = NULL, ...,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           environment = parent.frame()) {
    if (is.null(mapping)) {
    if (unit.out == "energy") {
      data <- q2e(data)
      mapping <- aes_(~w.length, ~s.e.response)
    } else if (unit.out %in% c("photon", "quantum")) {
      data <- e2q(data)
      mapping <- aes_(~w.length, ~s.q.response)
    } else {
      stop("Invalid 'unit.out' argument value: '", unit.out, "'")
    }
  }
  rmDerivedSpct(data)
  ggplot(data = data, mapping = mapping, ...,
         environment = environment)
}
#' @rdname ggplot
#'
#' @param plot.qty character string one of "transmittance" or "absorbance".
#'
#' @export
#'
ggplot.filter_spct <-
  function(data, mapping = NULL, ...,
           plot.qty = getOption("photobiology.filter.qty", default = "transmittance"),
           environment = parent.frame()) {
    if (is.null(mapping)) {
    if (plot.qty == "transmittance") {
      data <- A2T(data)
      mapping <- aes_(~w.length, ~Tfr)
    } else if (plot.qty == "absorbance") {
      data <- T2A(data)
      mapping <- aes_(~w.length, ~A)
    } else {
      stop("Invalid 'plot.qty' argument value: '", plot.qty, "'")
    }
  }
  rmDerivedSpct(data)
  ggplot(data = data, mapping = mapping, ...,
         environment = environment)
}
#' @rdname ggplot
#'
#' @note plot.qty is ignored for reflectors.
#'
#' @export
#'
ggplot.reflector_spct <- function(data, mapping = NULL, ..., plot.qty = NULL,
                               environment = parent.frame()) {
  if (is.null(mapping)) {
    mapping <- aes_(~w.length, ~Rfr)
  }
  rmDerivedSpct(data)
  ggplot(data = data, mapping = mapping, ...,
         environment = environment)
}
#' @rdname ggplot
#'
#' @export
#'
ggplot.cps_spct <- function(data, mapping = NULL, ...,
                                  environment = parent.frame()) {
  if (is.null(mapping)) {
    mapping <- aes_(~w.length, ~cps)
  }
  rmDerivedSpct(data)
  ggplot(data = data, mapping = mapping, ...,
         environment = environment)
}
#' @rdname ggplot
#'
#' @export
#'
ggplot.raw_spct <- function(data, mapping = NULL, ...,
                            environment = parent.frame()) {
  if (is.null(mapping)) {
    mapping <- aes_(~w.length, ~counts)
  }
  rmDerivedSpct(data)
  ggplot(data = data, mapping = mapping, ...,
         environment = environment)
}
