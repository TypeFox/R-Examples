

# Constructors ------------------------------------------------------------

#' Spectral-object constructor
#'
#' These fucntions can be used to create spectral objects derived from
#' \code{generic_spct}. They take as arguments numeric vectors for the data
#' character scalars for attributes, and a logical flag.
#'
#' @param w.length numeric vector with wavelengths in nanometres
#' @param s.e.irrad numeric vector with spectral energy irradiance in [W m-2
#'   nm-1] or [J d-1 m-2 nm-1]
#' @param s.q.irrad numeric A vector with spectral photon irradiance in [mol s-1
#'   m-2 nm-1] or [mol d-1 m-2 nm-1].
#' @param time.unit character string indicating the time unit used for spectral
#'   irradiance or exposure ("second" , "day" or "exposure") or an object of
#'   class duration as defined in package lubridate.
#' @param bswf.used character A string indicating the BSWF used, if any, for
#'   spectral effective irradiance or exposure ("none" or the name of the BSWF).
#' @param comment character A string to be added as a comment attribute to the
#'   object created.
#' @param strict.range logical Flag indicating whether off-range values result
#'   in an error instead of a warning.
#' @param ... other arguments passed to \code{data_frame()}
#'
#' @return A object of class generic_spct or a class derived from it, depending
#'   on the function used. In other words an object of a class with the same
#'   name as the constructor function.
#'
#' @export
#'
#' @family creation of spectral objects functions
#'
#' @rdname source_spct
#'
#' @note The functions can be used to add only one spectral quantity to a
#'   spectral object. Some of the functions have different arguments, for the
#'   same quantity expressed in different units. An actual parameter can be
#'   supplied to only one of these formal parameters in a given call to any
#'   of these functions.
#'
source_spct <- function(w.length = NULL, s.e.irrad = NULL, s.q.irrad = NULL,
                        time.unit = c("second", "day", "exposure"),
                        bswf.used = c("none", "unknown"),
                        comment = NULL, strict.range = FALSE, ...) {
  if (length(w.length) == 0) {
    z <- dplyr::data_frame(w.length = numeric(), s.e.irrad = numeric(), ...)
  } else if (is.null(s.q.irrad) && (is.numeric(s.e.irrad))) {
    z <- dplyr::data_frame(w.length, s.e.irrad, ...)
  } else if (is.null(s.e.irrad) && (is.numeric(s.q.irrad))) {
    z <- dplyr::data_frame(w.length, s.q.irrad, ...)
  } else {
    warning("Only one of s.e.irrad or s.q.irrad should be different from NULL.")
    z <- dplyr::data_frame(w.length, ...)
  }
  if (!is.null(comment)) {
    comment(z) <- comment
  }
  setSourceSpct(z,
                time.unit = time.unit,
                bswf.used = bswf.used,
                strict.range = strict.range)
  z
}

#' @rdname source_spct
#'
#' @param counts numeric vector with raw counts expressed per scan
#' @param instr.desc a list
#' @param instr.settings a list
#'
#' @export
#'
raw_spct <- function(w.length = NULL, counts = NA_real_, comment = NULL,
                     instr.desc = NA, instr.settings = NA, ...) {
  if (length(w.length) == 0) {
    z <- dplyr::data_frame(w.length = numeric(), counts = numeric(), ...)
  } else {
    z <- dplyr::data_frame(w.length = w.length, counts = counts, ...)
  }
  if (!is.null(comment)) {
    comment(z) <- comment
  }
  setRawSpct(z)
  setInstrDesc(z, instr.desc)
  setInstrSettings(z, instr.settings)
  z
}

#' @rdname source_spct
#'
#' @param cps numeric vector with linearized raw counts expressed per second
#'
#' @export
#'
cps_spct <- function(w.length = NULL, cps = NA_real_, comment = NULL,
                     instr.desc = NA, instr.settings = NA, ...) {
  if (any(grepl("^cps", names(list(...)))) && is.na(cps)) {
    if (length(w.length) == 0) {
      z <- dplyr::data_frame(w.length = numeric(), ...)
    } else {
      z <- dplyr::data_frame(w.length = w.length, ...)
    }
  } else {
    if (length(w.length) == 0) {
      z <- dplyr::data_frame(w.length = numeric(), cps = numeric(), ...)
    } else {
      z <- dplyr::data_frame(w.length = w.length, cps = cps, ...)
    }
  }
  if (!is.null(comment)) {
    comment(z) <- comment
  }
  setCpsSpct(z)
  setInstrDesc(z, instr.desc)
  setInstrSettings(z, instr.settings)
  z
}

#' @rdname source_spct
#'
#' @export
#'
generic_spct <- function(w.length = NULL, comment = NULL, ...) {
  if (length(w.length) == 0) {
    z <- dplyr::data_frame(w.length = numeric(), ...)
  } else {
    z <- dplyr::data_frame(w.length = w.length, ...)
  }
  if (!is.null(comment)) {
    comment(z) <- comment
  }
  setGenericSpct(z)
  z
}

#' @rdname source_spct
#'
#' @param s.e.response numeric vector with spectral energy irradiance in W m-2
#'   nm-1 or J d-1 m-2 nm-1
#' @param s.q.response numeric vector with spectral photon irradiance in mol s-1
#'   m-2 nm-1 or mol d-1 m-2 nm-1
#'
#' @export
#'
response_spct <- function(w.length = NULL, s.e.response = NULL, s.q.response = NULL,
                          time.unit = c("second", "day", "exposure"),
                          comment = NULL, ...) {
  if (length(w.length) == 0) {
    z <- dplyr::data_frame(w.length = numeric(), s.e.response = numeric(), ...)
  } else if (is.null(s.q.response) && (is.numeric(s.e.response))) {
    z <- dplyr::data_frame(w.length, s.e.response, ...)
  } else if (is.null(s.e.response) && (is.numeric(s.q.response))) {
    z <- dplyr::data_frame(w.length, s.q.response, ...)
  } else {
    warning("Only one of s.e.response or s.q.response should be different from NULL.")
    z <- dplyr::data_frame(w.length, ...)
  }
  if (!is.null(comment)) {
    comment(z) <- comment
  }
  setResponseSpct(z, time.unit)
  z
}

#' @rdname source_spct
#'
#' @param Tfr numeric vector with spectral transmittance as fraction of one
#' @param Tpc numeric vector with spectral transmittance as percent values
#' @param A   numeric vector of absorbance values (log10 based)
#' @param Tfr.type character string indicating whether transmittance values are
#'   "total" or "internal" values
#'
#' @note "internal" transmittance is defined as the transmittance of the
#'   material body itself, while "total" transmittance includes the effects of
#'   surface reflectance on the ammount of light transmitted.
#'
#' @export
#'
filter_spct <- function(w.length=NULL, Tfr=NULL, Tpc=NULL, A=NULL,
                        Tfr.type=c("total", "internal"),
                        comment=NULL, strict.range = FALSE, ...) {
  if (length(w.length) == 0) {
    z <- dplyr::data_frame(w.length = numeric(), Tfr = numeric())
  } else if (is.null(Tpc) && is.null(A) && is.numeric(Tfr)) {
    z <- dplyr::data_frame(w.length, Tfr, ...)
  } else if (is.null(Tfr) && is.null(A) && is.numeric(Tpc)) {
    z <- dplyr::data_frame(w.length, Tpc, ...)
  } else if (is.null(Tpc) && is.null(Tfr) && is.numeric(A)) {
    z <- dplyr::data_frame(w.length, A, ...)
  } else {
    warning("Only one of Tfr, Tpc or A should be different from NULL.")
    z <- dplyr::data_frame(w.length, ...)
  }
  if (!is.null(comment)) {
    comment(z) <- comment
  }
  setFilterSpct(z, Tfr.type, strict.range = strict.range)
  z
}

#' @rdname source_spct
#'
#' @param Rfr numeric vector with spectral refletance as fraction of one
#' @param Rpc numeric vector with spectral reflectance as percent values
#' @param Rfr.type character A string, either "total" or "specular".
#'
#' @export
#'
reflector_spct <- function(w.length = NULL, Rfr=NULL, Rpc=NULL,
                           Rfr.type=c("total", "specular"),
                           comment=NULL, strict.range = FALSE, ...) {
  if (length(w.length) == 0) {
    z <- dplyr::data_frame(w.length = numeric(), Rfr = numeric(), ...)
  } else if (is.null(Rpc) && is.numeric(Rfr)) {
    z <- dplyr::data_frame(w.length, Rfr, ...)
  } else if (is.null(Rfr) && is.numeric(Rpc)) {
    z <- dplyr::data_frame(w.length, Rpc, ...)
  } else {
    warning("Only one of Rfr, or Rpc should be different from NULL.")
    z <- dplyr::data_frame(w.length, ...)
  }
  if (!is.null(comment)) {
    comment(z) <- comment
  }
  setReflectorSpct(z, Rfr.type = Rfr.type, strict.range = strict.range)
  z
}

#' @rdname source_spct
#'
#' @export
#'
object_spct <- function(w.length=NULL, Rfr=NULL, Tfr=NULL,
                        Tfr.type=c("total", "internal"),
                        Rfr.type=c("total", "specular"),
                        comment=NULL, strict.range = FALSE, ...) {
  if (length(w.length) == 0) {
    z <- dplyr::data_frame(w.length = numeric(),
                           Rfr = numeric(), Tfr = numeric(), ...)
  } else {
    z <- dplyr::data_frame(w.length, Rfr, Tfr, ...)
  }
  if (!is.null(comment)) {
    comment(z) <- comment
  }
  setObjectSpct(z,
                Tfr.type = Tfr.type,
                Rfr.type = Rfr.type,
                strict.range = strict.range)
  z
}

#' @rdname source_spct
#'
#' @param x,y,z numeric colour coordinates
#'
#' @export
#'
chroma_spct <- function(w.length=NULL,
                        x,
                        y,
                        z,
                        comment=NULL, strict.range = FALSE, ...) {
  if (length(w.length) == 0) {
    z <- dplyr::data_frame(w.length = numeric(),
                           x = numeric(), y = numeric(), z = numeric(), ...)
  } else {
    z <- dplyr::data_frame(w.length, x, y, z, ...)
  if (!is.null(comment)) {
    comment(z) <- comment
  }
  setChromaSpct(z)
  }
  z
}

# as functions for spct classes --------------------------------------------

#' Spectral-object copy constructor
#'
#' Return a copy of an R object with its class set to a given type of spectrum.
#'
#' @param x an R object
#'
#' @return These functions return a copy of \code{x} converted into a given
#'   class of spectral object, if \code{x} is a valid argument to the
#'   correcponding set function.
#'
#' @export
#'
#' @family creation of spectral objects functions
#' @rdname as.generic_spct
#'
as.generic_spct <- function(x) {
  y <- x
  setGenericSpct(y)
}

#' @rdname as.generic_spct
#'
#' @export
#'
as.raw_spct <- function(x) {
  y <- x
  setRawSpct(y)
}

#' @rdname as.generic_spct
#'
#' @export
#'
as.cps_spct <- function(x) {
  y <- x
  setCpsSpct(y)
}

#' @rdname as.generic_spct
#'
#' @param time.unit character A string, "second", "day" or "exposure"
#' @param bswf.used character
#' @param strict.range logical Flag indicating whether off-range values result
#'   in an error instead of a warning
#'
#' @export
#'
as.source_spct <- function(x,
                           time.unit=c("second", "day", "exposure"),
                           bswf.used=c("none", "unknown"),
                           strict.range = FALSE) {
  y <- x
  setSourceSpct(y, time.unit, strict.range = strict.range, bswf.used = bswf.used)
}

#' @rdname as.generic_spct
#'
#' @export
#'
as.response_spct <- function(x, time.unit = "second") {
  y <- x
  setResponseSpct(y, time.unit = time.unit)
}

#' @rdname as.generic_spct
#'
#' @param Tfr.type a character string, either "total" or "internal"
#'
#' @export
#'
as.filter_spct <- function(x, Tfr.type=c("total", "internal"), strict.range = FALSE) {
  y <- x
  setFilterSpct(y, Tfr.type, strict.range = strict.range)
}

#' @rdname as.generic_spct
#'
#' @param Rfr.type a character string, either "total" or "specular"
#'
#' @export
#'
as.reflector_spct <- function(x, Rfr.type = c("total", "specular"), strict.range = FALSE) {
  y <- x
  setReflectorSpct(y, Rfr.type = Rfr.type, strict.range = strict.range)
}

#' @rdname as.generic_spct
#'
#' @export
#'
as.object_spct <- function(x,
                           Tfr.type=c("total", "internal"),
                           Rfr.type=c("total", "specular"),
                           strict.range = FALSE) {
  y <- x
  setObjectSpct(y, Tfr.type = Tfr.type, Rfr.type = Rfr.type,
                strict.range = strict.range)
}

#' @rdname as.generic_spct
#'
#' @export
#'
as.chroma_spct <- function(x) {
  y <- x
  setChromaSpct(y)
}


# merge -------------------------------------------------------------------


#' Merge two generic_spct objects
#'
#' Relatively quick merge of two spct objects based on w.length.
#'
#' @param x generic_spct (or derived) objects to be merged
#' @param y generic_spct (or derived) objects to be merged
#' @param by a vector of shared column names in \code{x} and \code{y} to merge on;
#' \code{by} defaults to \code{w.length}.
#' @param ... other arguments passed to \code{dplyr::inner_join()}
#'
#' @note if the class of x and y is the same, it is preserved, but
#' if it differs \code{generic_spct} is used for the returned value,
#' except when x and y, are one each of classes reflector_spct and
#' filter_spct in which case an object_spct is returned.
#' In the current implementation only wavelengths values shared
#' by x and y are preserved.
#'
#' @seealso \code{\link[dplyr]{join}}
#'
#' @export
#'
merge.generic_spct <- function(x, y, by = "w.length", ...) {
  class.x <- class(x)
  class.y <- class(y)
  if (identical(class.x, class.y)) {
    z <- dplyr::inner_join(x, y, by = by, ...)
    class(z) <- class.x
    warning("Attributes lost when merging two objects of class '", class.x, "'.")
  } else if ("filter_spct" %in% class.x && "reflector_spct" %in% class.y) {
    xx <- A2T(x, action = "replace", byref = FALSE)
    z <- dplyr::inner_join(xx, y, by = "w.length", ...)
    setObjectSpct(z, Tfr.type = getTfrType(x), Rfr.type = getRfrType(y))
  } else if ("reflector_spct" %in% class.x && "filter_spct" %in% class.y) {
    yy <- A2T(y, action = "replace", byref = FALSE)
    z <- dplyr::inner_join(xx, yy, by = "w.length", ...)
    setObjectSpct(z, Tfr.type = getTfrType(y), Rfr.type = getRfrType(x))
  } else {
    z <- dplyr::inner_join(x, y, by = "w.length", ...)
    setGenericSpct(z)
  }
  comment(z) <- paste("Merged spectrum\ncomment(x):\n",
                      comment(x),
                      "\nclass: ",
                      class_spct(x),
                      "\n\ncomment(y):\n",
                      comment(y),
                      "\nclass: ",
                      class_spct(y))
  z
}

