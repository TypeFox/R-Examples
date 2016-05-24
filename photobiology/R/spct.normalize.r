
# normalize ---------------------------------------------------------------

#' Normalize spectral data
#'
#' These functions return a spectral object of the same class as the one
#' supplied as argument but with the spectral data normalized to 1.o a certain
#' wavelength.
#'
#' @param x An R object
#' @param ... not used in current version
#' @return A new object of the same class as \code{x}.
#' @export normalize
#' @note Accepted values for \code{norm} vary depending on the class of
#'   \code{x}
#' @family rescaling functions
#'
normalize <- function(x, ...) UseMethod("normalize")

#' @describeIn normalize Default for generic function
#'
#' @export
normalize.default <- function(x, ...) {
  warning("'normalize' is not defined for objects of class ", class(x)[1])
  x
}


#' @describeIn normalize Normalize a \code{source_spct} object.
#'
#' @param range An R object on which \code{range()} returns a numeric vector of
#'   length 2 with the limits of a range of wavelengths in nm, with min annd max
#'   wavelengths (nm)
#' @param norm numeric Normalization wavelength (nm) or character string "max",
#'   or "min" for normalization at the corresponding wavelngth, or "integral" or
#'   "mean" for rescaling by dividing by these values.
#' @param unit.out character Allowed values "energy", and "photon",
#'   or its alias "quantum"
#'
#' @export
#'
normalize.source_spct <- function(x,
                                  ...,
                                  range = NULL,
                                  norm = "max",
                                  unit.out = getOption("photobiology.radiation.unit",
                                                       default="energy")) {
  if (unit.out == "energy") {
    return(normalize_spct(spct = q2e(x, action = "replace"),
                          range = range,
                          norm = norm,
                          col.names = "s.e.irrad"))
  } else if (unit.out %in% c("photon", "quantum") ) {
    return(normalize_spct(spct = e2q(x, action = "replace"),
                          range = range,
                          norm = norm,
                          col.names = "s.q.irrad"))
  } else {
    stop("'unit.out ", unit.out, " is unknown")
  }
}

#' @describeIn normalize Normalize a response spectrum.
#'
#' @export
#'
normalize.response_spct <- function(x,
                                  ...,
                                  range = NULL,
                                  norm = "max",
                                  unit.out = getOption("photobiology.radiation.unit",
                                                       default="energy")) {
  if (unit.out == "energy") {
    return(normalize_spct(spct = q2e(x, action = "replace"),
                          range = range,
                          norm = norm,
                          col.names = "s.e.response"))
  } else if (unit.out %in% c("photon", "quantum") ) {
    return(normalize_spct(spct = e2q(x, action = "replace"),
                          range = range,
                          norm = norm,
                          col.names = "s.q.response"))
  } else {
    stop("'unit.out ", unit.out, " is unknown")
  }
}

#' @describeIn normalize Normalize a filter spectrum.
#'
#' @param qty.out character string  Allowed values are "transmittance", and
#'   "absorbance" indicating on which quantity to apply the normalization.
#'
#' @export
#'
normalize.filter_spct <-
  function(x,
           ...,
           range = NULL,
           norm = "max",
           qty.out = getOption("photobiology.filter.qty",
                               default = "transmittance")) {
    if (qty.out == "transmittance") {
      return(normalize_spct(spct = A2T(x, action = "replace"),
                            range = range,
                            norm = norm,
                            col.names = "Tfr"))
    } else if (qty.out == "absorbance") {
      return(normalize_spct(spct = T2A(x, action = "replace"),
                            range = range,
                            norm = norm,
                            col.names = "A"))
    } else {
      stop("'qty.out ", qty.out, " is unknown")
    }
  }

#' @describeIn normalize Normalize a reflector spectrum.
#'
#' @export
#'
normalize.reflector_spct <-
  function(x,
           ...,
           range = NULL,
           norm = "max",
           qty.out = NULL) {
    normalize_spct(spct = x,
                   range = range,
                   norm = norm,
                   col.names = "Rfr")
  }

#' @describeIn normalize Normalize a raw spectrum.
#'
#' @export
#'
normalize.raw_spct <-
  function(x,
           ...,
           range = NULL,
           norm = "max") {
    normalize_spct(spct = x,
                   range = range,
                   norm = norm,
                   col.names = grep("^counts", names(x), value = TRUE))
  }

#' @describeIn normalize Normalize a cps spectrum.
#'
#' @export
#'
normalize.cps_spct <-
  function(x,
           ...,
           range = NULL,
           norm = "max") {
    normalize_spct(spct = x,
                   range = range,
                   norm = norm,
                   col.names = grep("^cps", names(x), value = TRUE))
  }

#' @describeIn normalize Normalize a raw spectrum.
#'
#' @param col.names character vector containing the names of columns or
#'   variables to which to apply the normalization.
#'
#' @export
#'
normalize.generic_spct <-
  function(x,
           ...,
           range = NULL,
           norm = "max",
           col.names) {

    normalize_spct(spct = x,
                   range = range,
                   norm = norm,
                   col.names = col.names)
  }

# collections of spectra --------------------------------------------------


#' @describeIn normalize Normalize the members of a source_mspct object.
#'
#' @export
#'
normalize.source_mspct <-
  function(x,
           ...,
           range = NULL,
           norm = "max",
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy")) {
    msmsply(x,
            normalize,
            range = range,
            norm = norm,
            unit.out = unit.out,
            ...)

  }

#' @describeIn normalize Normalize the members of a response_mspct object.
#'
#' @export
#'
normalize.response_mspct <-
  function(x,
           ...,
           range = NULL,
           norm = "max",
           unit.out = getOption("photobiology.radiation.unit",
                                default = "energy")) {
    msmsply(x,
            normalize,
            range = range,
            norm = norm,
            unit.out = unit.out,
            ...)

  }

#' @describeIn normalize Normalize the members of a filter_mspct object.
#'
#' @export
#'
normalize.filter_mspct <-
  function(x,
           ...,
           range = NULL,
           norm = "max",
           qty.out = getOption("photobiology.filter.qty",
                               default = "transmittance")) {
    msmsply(x,
            normalize,
            range = range,
            norm = norm,
            qty.out = qty.out,
            ...)

  }

#' @describeIn normalize Normalize the members of a reflector_mspct object.
#'
#' @export
#'
normalize.reflector_mspct <- function(x,
                                      ...,
                                      range = x,
                                      norm = "max",
                                      qty.out = NULL) {
  msmsply(x,
          normalize,
          range = range,
          norm = norm,
          qty.out = qty.out,
          ...)

}

#' @describeIn normalize Normalize the members of a raw_mspct object.
#'
#' @export
#'
normalize.raw_mspct <- function(x,
                                ...,
                                range = x,
                                norm = "max") {
  msmsply(x,
          normalize,
          range = range,
          norm = norm,
          ...)
}

#' @describeIn normalize Normalize the members of a cps_mspct object.
#'
#' @export
#'
normalize.cps_mspct <- function(x,
                                ...,
                                range = x,
                                norm = "max") {
  msmsply(x,
          normalize,
          range = range,
          norm = norm,
          ...)
}

# PRIVATE -----------------------------------------------------------------

#' @keywords internal
#'
normalize_spct <- function(spct, range, norm, col.names) {
  stopifnot(is.any_spct(spct), !is.null(col.names),
            col.names %in% names(spct))
  if (is.null(range) || all(is.na(range))) {
    range <- range(spct)
  }
  tmp.spct <- trim_spct(spct, range)
  # rescaling needed
  if (!is.null(norm)) {
    for (col in col.names) {

      if (is.character(norm)) {
        if (norm %in% c("max", "maximum")) {
          idx <- which.max(tmp.spct[[col]])
        } else if (norm %in% c("min", "minimum")) {
          idx <- which.min(tmp.spct[[col]])
        } else {
          warning("Invalid character '", norm, "'value in 'norm'")
          idx <- NA
        }
        scale.factor <- 1 / tmp.spct[idx, col, drop = TRUE]
        norm <- tmp.spct[idx, "w.length", drop = TRUE]
      } else if (is.numeric(norm)) {
        if (norm >= min(tmp.spct) && norm <= max(tmp.spct)) {
          tmp.spct <- tmp.spct[ , c("w.length", col)]
          class(tmp.spct) <- class(spct)
          scale.factor <- 1 /
            interpolate_spct(spct = tmp.spct, w.length.out = norm)[ , eval(col)]
        } else {
          warning("'norm = ", norm, "' value outside spectral data range of ",
                  round(min(tmp.spct), 1), " to ", round(max(tmp.spct), 1), " (nm)")
          scale.factor <- NA
        }
      } else {
        stop("'norm' should be numeric or character")
      }
      spct[[col]] <- spct[ , col, drop = TRUE] * scale.factor
    }
  } else {
    return(spct)
  }
  spct <- setNormalized(spct, norm)
  spct
}


# is_normalized function --------------------------------------------------

#' Query whether a generic spectrum has been normalized.
#'
#' This function tests a \code{generic_spct} object for an attribute that
#' signals whether the spectral data has been normalized or not after the object
#' was created.
#'
#' @param x An R object.
#'
#' @return A \code{logical} value. If \code{x} is not normalized or \code{x} is
#'   not a \code{generic_spct} object the value returned is \code{FALSE}.
#'
#' @export
#' @family rescaling functions
#'
is_normalized <- function(x) {
  if (!is.any_spct(x) && !is.any_summary_spct(x)) {
    return(NA)
  }
  spct.attr <- attr(x, "normalized", exact = TRUE)
  as.logical(!is.null(spct.attr) && as.logical(spct.attr))
}

# getNormalized -----------------------------------------------------------

#' Get the "normalized" attribute
#'
#' Funtion to read the "normalized" attribute of an existing generic_spct
#' object.
#'
#' @param x a generic_spct object
#'
#' @return character or numeric or logical
#'
#' @note if x is not a \code{generic_spct} object, \code{NA} is returned
#'
#' @export
#' @family rescaling functions
#'
getNormalized <- function(x) {
  if (is.any_spct(x) || is.any_summary_spct(x)) {
    normalized <- attr(x, "normalized", exact = TRUE)
    if (is.null(normalized) || is.na(normalized)) {
      # need to handle objects created with old versions
      normalized <- FALSE
    }
    return(normalized[[1]])
  } else {
    return(NA)
  }
}

#' Set the "normalized" attribute
#'
#' Funtion to write the "normalized" attribute of an existing generic_spct
#' object.
#'
#' @param x a generic_spct object
#' @param norm numeric or logical
#'
#' @note if x is not a \code{generic_spct} object, x is not modified.
#'
#' @export
#' @family rescaling functions
#'
setNormalized <- function(x, norm = FALSE) {
  name <- substitute(x)
  if ((is.any_spct(x) || is.any_summary_spct(x)) &&
      (is.na(norm) || is.numeric(norm) || is.logical(norm))) {
    attr(x, "normalized") <- norm
    if (is.name(name)) {
      name <- as.character(name)
      assign(name, x, parent.frame(), inherits = TRUE)
    }
  }
  invisible(x)
}

