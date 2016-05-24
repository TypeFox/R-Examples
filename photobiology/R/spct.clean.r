#' Clean (=replace) off-range values in a spectrum
#'
#' These functions implement the equivalent of replace() but for spectral
#' objects instead of vectors.
#'
#' @param x an R object
#' @param range numeric vector of wavelengths
#' @param range.s.data numeric vector of length two giving the allowable
#'     range for the spectral data.
#' @param fill numeric vector of length 1 or 2, giving the replacement
#'     values to use at each extreme of the range.
#' @param ... currently ignored
#'
#'
#' @export
#'
clean <- function(x, range, range.s.data, fill, ...) UseMethod("clean")

#' @describeIn clean Default for generic function
#'
#' @export
#'
clean.default <- function(x, range, range.s.data, fill, ...) {
  warning("'clean()' is not defined for objects of class '", class(x)[1], "'.")
  x
}

#' @describeIn clean Replace off-range values in a source spectrum
#' @param unit.out character string with allowed values "energy", and "photon",
#'   or its alias "quantum"
#'
#' @export
#'
clean.source_spct <-
  function(x,
           range = x,
           range.s.data = c(0,NA),
           fill = range.s.data,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...) {
    if (unit.out == "quantum" || unit.out == "photon") {
      clean_spct(x = e2q(x, action = "replace"),
                 range = range,
                 range.s.data = range.s.data,
                 fill = fill,
                 col.names = "s.q.irrad",
                 ...)
    } else if (unit.out == "energy") {
      clean_spct(x = q2e(x, action = "replace"),
                 range = range,
                 range.s.data = range.s.data,
                 fill = fill,
                 col.names = "s.e.irrad",
                 ...)
    } else {
      stop("unit.out: '", unit.out, "' not supported")
    }
  }

#' @describeIn clean Replace off-range values in a filter spectrum
#' @param qty.out character string with allowed values "energy", and "photon",
#'   or its alias "quantum"
#'
#' @export
#'
clean.filter_spct <-
  function(x,
           range = x,
           range.s.data = NULL,
           fill = range.s.data,
           qty.out = getOption("photobiology.filter.qty", default = "transmittance"),
           ...) {
    if (is.null(range.s.data)) {
      if (qty.out == "transmittance") {
        range.s.data <- c(0,1)
      } else {
        range.s.data <- c(0,NA)
      }
      fill <- range.s.data
    }
    if (qty.out == "transmittance") {
      clean_spct(x = A2T(x, action = "replace"),
                 range = range,
                 range.s.data = range.s.data,
                 fill = fill,
                 col.names = "Tfr",
                 ...)
    } else if (qty.out == "absorbance") {
      clean_spct(x = T2A(x, action = "replace"),
                 range = range,
                 range.s.data = range.s.data,
                 fill = fill,
                 col.names = "A",
                 ...)
    } else {
      stop("qty.out: '", qty.out, "' not supported")
    }
  }

#' @describeIn clean Replace off-range values in a reflector spectrum
#'
#' @export
#'
clean.reflector_spct <-
  function(x,
           range = x,
           range.s.data = c(0, 1),
           fill = range.s.data,
           ...) {
    clean_spct(x = x,
               range = range,
               range.s.data = range.s.data,
               fill = fill,
               col.names = "Rfr",
               ...)
  }


#' @describeIn clean Replace off-range values in a response spectrum
#'
#' @export
#'
clean.response_spct <-
  function(x,
           range = x,
           range.s.data = c(0,NA),
           fill = range.s.data,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           ...) {
    if (unit.out == "quantum" || unit.out == "photon") {
      clean_spct(x = e2q(x, action = "replace"),
                 range = range,
                 range.s.data = range.s.data,
                 fill = fill,
                 col.names = "s.q.response",
                 ...)
    } else if (unit.out == "energy") {
      clean_spct(x = q2e(x, action = "replace"),
                 range = range,
                 range.s.data = range.s.data,
                 fill = fill,
                 col.names = "s.e.response",
                 ...)
    } else {
      stop("unit.out: '", unit.out, "' not supported")
    }
  }

#' @describeIn clean Replace off-range values in a counts per second spectrum
#'
#' @export
#'
clean.cps_spct <-
  function(x,
           range = x,
           range.s.data = c(0, NA),
           fill = range.s.data,
           ...) {
    clean_spct(x = x,
               range = range,
               range.s.data = range.s.data,
               fill = fill,
               col.names = grep("^cps", names(x), value = TRUE),
               ...)
  }

#' @describeIn clean Replace off-range values in a raw counts spectrum
#'
#' @export
#'
clean.raw_spct <-
  function(x,
           range = x,
           range.s.data = c(NA_real_, NA_real_),
           fill = range.s.data,
           ...) {
    clean_spct(x = x,
               range = range,
               range.s.data = range.s.data,
               fill = fill,
               col.names = grep("^counts", names(x), value = TRUE),
               ...)
  }

#' @describeIn clean Replace off-range values in a generic spectrum
#'
#' @param col.names character The name of the variable to clean
#'
#' @export
#'
clean.generic_spct <-
  function(x,
           range = x,
           range.s.data = c(NA_real_, NA_real_),
           fill = range.s.data,
           col.names,
           ...) {
    clean_spct(x = x,
               range = range,
               range.s.data = range.s.data,
               fill = fill,
               col.names = col.names,
               ...)
  }

# Collections of spectra --------------------------------------------------

#' @describeIn clean
#'
#' @export
#'
clean.source_mspct <-
  function(x,
           range = NULL,
           range.s.data = c(0,NA),
           fill = range.s.data,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...) {
    if (is.null(range)) {
      msmsply(mspct = x,
              .fun = clean,
              range.s.data = range.s.data,
              fill = fill,
              unit.out = unit.out,
              ...)
    } else {
      msmsply(mspct = x,
              .fun = clean,
              range = range,
              range.s.data = range.s.data,
              fill = fill,
              unit.out = unit.out,
              ...)
    }
  }

#' @describeIn clean
#'
#' @export
#'
clean.filter_mspct <-
  function(x,
           range = NULL,
           range.s.data = NULL,
           fill = range.s.data,
           qty.out = getOption("photobiology.filter.qty",
                               default = "transmittance"),
           ...) {
    if (is.null(range)) {
      msmsply(mspct = x,
              .fun = clean,
              range.s.data = range.s.data,
              fill = fill,
              qty.out = qty.out,
              ...)
    } else {
      msmsply(mspct = x,
              .fun = clean,
              range = range,
              range.s.data = range.s.data,
              fill = fill,
              qty.out = qty.out,
              ...)
    }
   }

#' @describeIn clean
#'
#' @export
#'
clean.reflector_mspct <-
  function(x,
           range = NULL,
           range.s.data = c(0, 1),
           fill = range.s.data,
           ...) {
    if (is.null(range)) {
      msmsply(mspct = x,
              .fun = clean,
              range.s.data = range.s.data,
              fill = fill,
              ...)
    } else {
      msmsply(mspct = x,
              .fun = clean,
              range = range,
              range.s.data = range.s.data,
              fill = fill,
              ...)
    }
   }

#' @describeIn clean
#'
#' @export
#'
clean.response_mspct <-
  function(x,
           range = NULL,
           range.s.data = c(0,NA),
           fill = range.s.data,
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy"),
           ...) {
    if (is.null(range)) {
      msmsply(mspct = x,
              .fun = clean,
              range.s.data = range.s.data,
              fill = fill,
              unit.out = unit.out,
              ...)
    } else {
      msmsply(mspct = x,
              .fun = clean,
              range = range,
              range.s.data = range.s.data,
              fill = fill,
              unit.out = unit.out,
              ...)
    }
  }

#' @describeIn clean
#'
#' @export
#'
clean.cps_mspct <-
  function(x,
           range = NULL,
           range.s.data = c(0, NA),
           fill = range.s.data,
           ...) {
    if (is.null(range)) {
      msmsply(mspct = x,
              .fun = clean,
              range.s.data = range.s.data,
              fill = fill,
              ...)
    } else {
      msmsply(mspct = x,
              .fun = clean,
              range = range,
              range.s.data = range.s.data,
              fill = fill,
              ...)
    }
  }

#' @describeIn clean
#'
#' @export
#'
clean.raw_mspct <- clean.cps_mspct

#' @describeIn clean
#'
#' @export
#'
clean.generic_mspct <-
  function(x,
           range = x,
           range.s.data = c(NA_real_, NA_real_),
           fill = range.s.data,
           col.names,
           ...) {
    if (is.null(range)) {
      msmsply(mspct = x,
              .fun = clean,
              range.s.data = range.s.data,
              fill = fill,
              ...)
    } else {
      msmsply(mspct = x,
              .fun = clean,
              range = range,
              range.s.data = range.s.data,
              fill = fill,
              ...)
    }
  }

# PRIVATE -----------------------------------------------------------------

#' Clean a spectrum
#'
#' These functions implement the equivalent of replace() but for spectral
#' objects instead of vectors.
#'
#' @param x an R object
#' @param range numeric vector of wavelengths
#' @param range.s.data numeric vector of length two giving the allowable
#'     range for the spectral data.
#' @param fill numeric vector of length 1 or 2, giving the replacement
#'     values to use at each extreme of the range.
#' @param col.names character The name of the variable to clean.
#' @param ... currently ignored
#'
#' @keywords internal
#'
clean_spct <-
  function(x,
           range,
           range.s.data,
           fill,
           col.names,
           col.pattern = NULL,
           ...) {
    stopifnot(length(range) >= 2L &&
                length(range.s.data) == 2L &&
                length(fill) <= 2L)
    if (length(fill) == 1) {
      fill <- c(fill, fill)
    }
    # wavelength range
    if (is.any_spct(range) || is.numeric(range) && length(range) > 2L) {
      range <- range(range, na.rm = TRUE)
    } else {
      if (is.na(range[1])) {
        range[1] <- min(x)
      }
      if (is.na(range[2])) {
        range[2] <- max(x)
      }
    }
    selector <- x[["w.length"]] >= range[1] & x[["w.length"]] <= range[2]

    if (is.na(range.s.data[1])) {
      range.s.data[1] <- -Inf
    }
    if (is.na(range.s.data[2])) {
      range.s.data[2] <- Inf
    }

    for (col in col.names) {
      x[selector, col] <-  ifelse(x[selector, col] < range.s.data[1],
                                  fill[1],
                                  ifelse(x[selector, col] > range.s.data[2],
                                         fill[2],
                                         x[selector, col]))
    }
    x
  }

