#' Find peaks in a spectrum
#'
#' This function finds all peaks (local maxima) in a spectrum, using a user
#' selectable size threshold relative to the tallest peak (global maximum). This
#' a wrapper built on top of function peaks from package splus2R.
#'
#' @param x numeric array
#' @param ignore_threshold numeric value between 0.0 and 1.0 indicating the size
#'   threshold below which peaks will be ignored.
#' @param span a peak is defined as an element in a sequence which is greater
#'   than all other elements within a window of width span centered at that
#'   element. The default value is 3, meaning that a peak is bigger than both of
#'   its neighbors. Default: 3.
#' @param strict ogical flag: if TRUE, an element must be strictly greater than
#'   all other values in its window to be considered a peak. Default: TRUE.
#'
#' @return an object like s.irrad of logical values. Values that are TRUE
#'   correspond to local peaks in the data.
#'
#' @export
#' @examples
#' with(sun.data, w.length[find_peaks(s.e.irrad)])
#'
#' @note This function is a wrapper built on function
#'   \code{\link[splus2R]{peaks}} from \pkg{splus2R} and handles non-finite
#'   (including NA) values differently than \code{peaks}, instead of giving an
#'   error they are replaced with the smallest finite value in \code{x}.
#'
#' @seealso \code{\link[splus2R]{peaks}}
#'
#' @family peaks and valleys functions
#'
find_peaks <-
  function(x,
           ignore_threshold = 0.0,
           span = 3,
           strict = TRUE) {
    range_x <- range(x, finite = TRUE)
    min_x <- range_x[1]
    max_x <- range_x[2]
    x <- ifelse(!is.finite(x), min_x, x)
    # the next two lines catter for the case when max_x < 0, which is quite common with logs
    delta <- max_x - min_x
    top_flag <- ignore_threshold > 0.0
    scaled_threshold <- delta * abs(ignore_threshold)
    pks <- splus2R::peaks(x = x, span = span, strict = strict)
    if (abs(ignore_threshold) < 1e-5)
      return(pks)
    if (top_flag) {
      return(ifelse(x - min_x > scaled_threshold, pks , FALSE))
    } else {
      return(ifelse(max_x - x > scaled_threshold, pks , FALSE))
    }
  }

#' Get peaks and valleys in a spectrum
#'
#' These functions find peaks (local maxima) or valleys (local minima) in a
#' spectrum, using a user selectable size threshold relative to the tallest peak
#' (global maximum). This a wrapper built on top of function peaks from package
#' splus2R.
#'
#' @param x numeric
#' @param y numeric
#' @param ignore_threshold numeric Value between 0.0 and 1.0 indicating the
#'   relative size compared to tallest peak or deepest valley of the peaks
#'   to return.
#' @param span numeric A peak is defined as an element in a sequence which is
#'   greater than all other elements within a window of width \code{span}
#'   centered at that element. For example, a value of 3 means that a peak is
#'   bigger than both of its neighbors.
#' @param strict logical Flag: if TRUE, an element must be strictly greater than
#'   all other values in its window to be considered a peak. Default: TRUE.
#' @param x_unit character Vector of texts to be pasted at end of labels built
#'   from x value at peaks.
#' @param x_digits numeric Number of significant digits in wavelength label.
#'
#' @return A data frame with variables w.length and s.irrad with their values at
#'   the peaks or valleys plus a character variable of labels.
#'
#' @export
#' @examples
#' with(sun.spct, get_peaks(w.length, s.e.irrad))
#' with(sun.spct, get_valleys(w.length, s.e.irrad))
#'
#' @family peaks and valleys functions
#'
get_peaks <- function(x,
                      y,
                      ignore_threshold = 0.0,
                      span = 5,
                      strict = TRUE,
                      x_unit = "",
                      x_digits = 3) {
  stopifnot(length(x) == length(y))
  selector <- find_peaks(y, ignore_threshold, span, strict)
  if (sum(selector) < 1) {
    return(data.frame(
      x = numeric(0),
      y = numeric(0),
      label = character(0)
    ))
  } else {
    peaks.x <- x[selector]
    peaks.y <- y[selector]
    return(data.frame(
      x = peaks.x,
      y = peaks.y,
      label = paste(as.character(signif(
        x = peaks.x, digits = x_digits
      )), x_unit, sep = "")
    ))
  }
}

#' @rdname get_peaks
#' @export
#'
get_valleys <- function(x, y,
                        ignore_threshold = 0.0,
                        span = 5,
                        strict = TRUE,
                        x_unit = "",
                        x_digits = 3) {
  xy.data <- get_peaks(x, -y,
                       -ignore_threshold,
                       span = span,
                       strict = strict,
                       x_unit = x_unit,
                       x_digits = x_digits)
  xy.data$y <- -xy.data$y
  return(xy.data)
}


# peaks -------------------------------------------------------------------

#' Peaks or local maxima
#'
#' Function that returns a subset of an R object with observations corresponding
#' to local maxima.
#'
#' @param x an R object
#' @param ignore_threshold numeric value between 0.0 and 1.0 indicating the
#'   relative size compared to talelst peakthreshold below which peaks will be
#'   ignored.
#' @param span a peak is defined as an element in a sequence which is greater
#'   than all other elements within a window of width span centered at that
#'   element. The default value is 3, meaning that a peak is bigger than both of
#'   its neighbors. Default: 3.
#' @param strict logical flag: if TRUE, an element must be strictly greater than
#'   all other values in its window to be considered a peak. Default: TRUE.
#' @param ... ignored
#'
#' @return a subset of x with rows corresponding to local maxima.
#'
#'
#' @export
#'
#' @family peaks and valleys functions
#'
peaks <- function(x, span, ignore_threshold, strict, ...) UseMethod("peaks")

#' @describeIn peaks Default function usable on numeric vectors.
#' @export
peaks.default <- function(x, span, ignore_threshold, strict, ...) {
  x[NA]
}

#' @describeIn peaks Default function usable on numeric vectors.
#' @export
peaks.numeric <- function(x, span = 5, ignore_threshold, strict = TRUE, ...) {
  splus2R::peaks(x = x, span = span, strict = strict)
}

#' @describeIn peaks  Method for "generic_spct" objects.
#'
#' @export
#'
peaks.generic_spct <- function(x, span, ignore_threshold, strict, ...) {
  peaks.idx <- find_peaks(x[[names(x)[2]]],
                          span = span, ignore_threshold = ignore_threshold,
                          strict = strict)
  x[peaks.idx, ]
}

#' @describeIn peaks  Method for "source_spct" objects.
#'
#' @param unit.out character One of "energy" or "photon"
#'
#' @export
#'
#' @examples
#' peaks(sun.spct)
#'
peaks.source_spct <-
  function(x, span = 5, ignore_threshold = 0.0, strict = TRUE,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           ...) {
    if (unit.out == "energy") {
      z <- q2e(x, "replace", FALSE)
      col.name <- "s.e.irrad"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, "replace", FALSE)
      col.name <- "s.q.irrad"
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }
    peaks.idx <- find_peaks(z[[col.name]],
                            span = span, ignore_threshold = ignore_threshold,
                            strict = strict)
    z[peaks.idx, ]
  }

#' @describeIn peaks  Method for "response_spct" objects.
#'
#' @export
#'
peaks.response_spct <-
  function(x, span = 5, ignore_threshold = 0.0, strict = TRUE,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           ...) {
    if (unit.out == "energy") {
      z <- q2e(x, "replace", FALSE)
      col.name <- "s.e.response"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, "replace", FALSE)
      col.name <- "s.q.response"
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }
    peaks.idx <- find_peaks(z[[col.name]],
                            span = span, ignore_threshold = ignore_threshold,
                            strict = strict)
    z[peaks.idx, ]
  }

#' @describeIn peaks  Method for "filter_spct" objects.
#'
#' @param filter.qty character One of "transmittance" or "absorbance"
#'
#' @export
#'
peaks.filter_spct <-
  function(x, span = 5, ignore_threshold = 0, strict = TRUE,
           filter.qty = getOption("photobiology.filter.qty", default = "transmittance"),
           ...) {
    if (filter.qty == "transmittance") {
      z <- A2T(x, "replace", FALSE)
      col.name <- "Tfr"
    } else if (filter.qty == "absorbance") {
      z <- T2A(x, "replace", FALSE)
      col.name <- "A"
    } else {
      stop("Unrecognized 'filter.qty': ", filter.qty)
    }
    peaks.idx <- find_peaks(z[[col.name]],
                            span = span, ignore_threshold = ignore_threshold,
                            strict = strict)
    z[peaks.idx, ]
  }

#' @describeIn peaks  Method for "reflector_spct" objects.
#'
#' @export
#'
peaks.reflector_spct <- function(x, span = 5, ignore_threshold = 0, strict = TRUE, ...) {
  peaks.idx <- find_peaks(x[["Rfr"]],
                          span = span, ignore_threshold = ignore_threshold,
                          strict = strict)
  subset(x, idx = peaks.idx)
}

#' @describeIn peaks  Method for "cps_spct" objects.
#'
#' @export
#'
peaks.cps_spct <- function(x, span = 5, ignore_threshold = 0, strict = TRUE, ...) {
  peaks.idx <- find_peaks(x[["cps"]],
                          span = span, ignore_threshold = ignore_threshold,
                          strict = strict)
  x[peaks.idx, ]
}

#' @describeIn peaks  Method for "cps_spct" objects.
#'
#' @export
#'
peaks.generic_mspct <- function(x, span = 5, ignore_threshold = 0, strict = TRUE, ...) {
  msmsply(x,
          .fun = peaks,
          span = span,
          ignore_threshold = ignore_threshold,
          strict = strict,
          ... )
  }

# valleys -------------------------------------------------------------------

#' Valleys or local minima
#'
#' Function that returns a subset of an R object with observations corresponding
#' to local maxima.
#'
#' @param x an R object
#' @param ignore_threshold numeric value between 0.0 and 1.0 indicating the
#'   relative size compared to talelst peakthreshold below which valleys will be
#'   ignored.
#' @param span a peak is defined as an element in a sequence which is greater
#'   than all other elements within a window of width span centered at that
#'   element. The default value is 3, meaning that a peak is bigger than both of
#'   its neighbors. Default: 3.
#' @param strict logical flag: if TRUE, an element must be strictly greater than
#'   all other values in its window to be considered a peak. Default: TRUE.
#' @param ... ignored
#'
#' @return a subset of x with rows corresponding to local maxima.
#'
#'
#' @export
#'
#' @family peaks and valleys functions
#'
valleys <- function(x, span, ignore_threshold, strict, ...) UseMethod("valleys")

#' @describeIn valleys Default function usable on numeric vectors.
#' @export
valleys.default <- function(x, span, ignore_threshold, strict, ...) {
  x[which.max(x)]
}

#' @describeIn valleys Default function usable on numeric vectors.
#' @export
valleys.numeric <- function(x, span = 5, ignore_threshold, strict = TRUE, ...) {
  x[splus2R::peaks(x = -x, span = span, strict = strict)]
}

#' @describeIn valleys  Method for "generic_spct" objects.
#'
#' @export
#'
valleys.generic_spct <- function(x, span = 5, ignore_threshold = 0.0, strict = TRUE, ...) {
  valleys.idx <- find_peaks(-x[names(x)[2]],
                          span = span, ignore_threshold = ignore_threshold,
                          strict = strict)
  x[valleys.idx, ]
}

#' @describeIn valleys  Method for "source_spct" objects.
#'
#' @param unit.out character One of "energy" or "photon"
#'
#' @export
#'
#' @examples
#' valleys(sun.spct)
#'
valleys.source_spct <-
  function(x, span = 5, ignore_threshold = 0.0, strict = TRUE,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           ...) {
    if (unit.out == "energy") {
      z <- q2e(x, "replace", FALSE)
      col.name <- "s.e.irrad"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, "replace", FALSE)
      col.name <- "s.q.irrad"
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }
    valleys.idx <- find_peaks(-z[[col.name]],
                          span = span, ignore_threshold = ignore_threshold,
                          strict = strict)
    z[valleys.idx, ]
  }

#' @describeIn valleys  Method for "response_spct" objects.
#'
#' @export
#'
valleys.response_spct <-
  function(x, span = 5, ignore_threshold = 0.0, strict = TRUE,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           ...) {
    if (unit.out == "energy") {
      z <- q2e(x, "replace", FALSE)
      col.name <- "s.e.response"
    } else if (unit.out %in% c("photon", "quantum")) {
      z <- e2q(x, "replace", FALSE)
      col.name <- "s.q.response"
    } else {
      stop("Unrecognized 'unit.out': ", unit.out)
    }
    valleys.idx <- find_peaks(-z[[col.name]],
                            span = span, ignore_threshold = ignore_threshold,
                            strict = strict)
    z[valleys.idx, ]
  }

#' @describeIn valleys  Method for "filter_spct" objects.
#'
#' @param filter.qty character One of "transmittance" or "absorbance"
#'
#' @export
#'
valleys.filter_spct <-
  function(x, span = 5, ignore_threshold = 0, strict = TRUE,
           filter.qty = getOption("photobiology.filter.qty", default = "transmittance"),
           ...) {
    if (filter.qty == "transmittance") {
      z <- A2T(x, "replace", FALSE)
      col.name <- "Tfr"
    } else if (filter.qty == "absorbance") {
      z <- T2A(x, "replace", FALSE)
      col.name <- "A"
    } else {
      stop("Unrecognized 'filter.qty': ", filter.qty)
    }
    valleys.idx <- find_peaks(-z[[col.name]],
                              span = span, ignore_threshold = ignore_threshold,
                              strict = strict)
    z[valleys.idx, ]
  }

#' @describeIn valleys  Method for "reflector_spct".
#'
#' @export
#'
valleys.reflector_spct <- function(x, span = 5, ignore_threshold = 0, strict = TRUE, ...) {
  valleys.idx <- find_peaks(-x[["Rfr"]],
                          span = span, ignore_threshold = ignore_threshold,
                          strict = strict)
  x[valleys.idx, ]
}

#' @describeIn valleys  Method for "cps_spct" objects.
#'
#' @export
#'
valleys.cps_spct <- function(x, span = 5, ignore_threshold = 0, strict = TRUE, ...) {
  valleys.idx <- find_peaks(-x[["cps"]],
                          span = span, ignore_threshold = ignore_threshold,
                          strict = strict)
  x[valleys.idx, ]
}

#' @describeIn valleys  Method for "generic_mspct" objects.
#'
#' @export
#'
valleys.generic_mspct <- function(x, span = 5, ignore_threshold = 0, strict = TRUE, ...) {
  msmsply(x,
          .fun = valleys,
          span = span,
          ignore_threshold = ignore_threshold,
          strict = strict,
          ... )
}
