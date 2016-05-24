#' Trim (or expand) head and/or tail
#'
#' Trimming of head and tail of a spectrum based on wavelength limits,
#' interpolating the values at the boundaries. Trimming is needed for example to
#' remove short wavelength noise when the measured spectrum extends beyond the
#' known emission spectrum of the measured light source. Occasionally one may
#' want also to expand the wavelength range.
#'
#' @param spct an object of class "generic_spct"
#' @param range a numeric vector of length two, or any other object for which
#'   function range() will return two
#' @param low.limit shortest wavelength to be kept (defaults to shortest
#'   w.length value)
#' @param high.limit longest wavelength to be kept (defaults to longest w.length
#'   value)
#' @param use.hinges logical, if TRUE (the default) wavelengths in nm.
#' @param fill if fill==NULL then tails are deleted, otherwise tails or s.irrad
#'   are filled with the value of fill
#' @param byref logical indicating if new object will be created by reference or
#'   by copy of spct
#' @param verbose logical
#'
#' @return a spectrum of same class as input with its tails trimmed or expanded
#'
#' @note When expanding an spectrum, if fill==NULL, then expansion is not
#'   performed. Range can be "waveband" object, a numeric vector or a list of
#'   numeric vectors, or any other user-defined or built-in object for which
#'   \code{range()} returns a numeric vector of legth two, that can be
#'   interpreted as wavelengths expressed in nm.
#' @family trim functions
#'
#' @export
#' @examples
#' trim_spct(sun.spct, low.limit=300)
#' trim_spct(sun.spct, low.limit=300, fill=NULL)
#' trim_spct(sun.spct, low.limit=300, fill=NA)
#' trim_spct(sun.spct, low.limit=300, fill=0.0)
#' trim_spct(sun.spct, range = c(300, 400))
#' trim_spct(sun.spct, range = c(300, NA))
#' trim_spct(sun.spct, range = c(NA, 400))
trim_spct <- function(spct,
                      range = NULL,
                      low.limit = NULL, high.limit = NULL,
                      use.hinges = TRUE,
                      fill = NULL,
                      byref = FALSE,
                      verbose = getOption("photobiology.verbose", default = TRUE) )
{
  if (nrow(spct) == 0) {
    return(spct)
  }
  stopifnot(is.any_spct(spct))
  x <- spct
  if (is.null(use.hinges)) {
    use.hinges <- auto_hinges(spct[["w.length"]])
  }
  if (byref) {
    name <- substitute(spct)
  }
  class_spct <- class(spct)
  if (is.null(range)) {
    range[1] <- ifelse(is.null(low.limit), NA, low.limit)
    range[2] <- ifelse(is.null(high.limit), NA, high.limit)
  }
  range <- normalize_range_arg(range)
  low.limit <- range[1]
  high.limit <- range[2]
  trim.low <- low.limit > min(spct, na.rm = TRUE)
  trim.high <- high.limit < max(spct, na.rm = TRUE)
  if (trim.low && trim.high && high.limit - low.limit < 1e-7) {
    warning("When trimming, 'range' must be a finite wavelength interval > 1E-7 nm")
    return(spct[FALSE, ]) # returns a spct object with nrow equal to zero
  }
  names.spct <- names(spct)
  names.data <- names.spct[names.spct != "w.length"]
  # comment.spct <- comment(spct)
  # time.unit.spct <- getTimeUnit(spct)
  # Tfr.type.spct <- getTfrType(spct)
  # Rfr.type.spct <- getRfrType(spct)
  # when.measured.spct <- getWhenMeasured(spct)
  # where.measured.spct <- getWhereMeasured(spct)
  # what.measured.spct <- getWhatMeasured(spct)
  # instr.settings.spct <- getInstrSettings(spct)
  # check whether we should expand the low end
  low.end <- min(spct, na.rm = TRUE)
  if (trim.low && low.end > low.limit) {
    if (!is.null(fill)) {
      # expand short tail
      low.tail.length <-  trunc(low.end - low.limit) + ifelse(use.hinges, 2, 1)
      low.tail.w.length <- seq(from = low.limit,
                               to = ifelse(use.hinges, low.end - 1e-12, low.end - 1),
                               length.out = low.tail.length)
      spct.top <- dplyr::data_frame(w.length = low.tail.w.length)
      for (data.col in names.data) {
        spct.top[[data.col]] <- fill
      }
      spct <- plyr::rbind.fill(list(spct.top, spct))
      spct <- dplyr::as_data_frame(spct)
      setGenericSpct(spct)
      low.end <- min(spct)
    } else {
        if (verbose && (low.end - low.limit) > 0.01) {
          warning("Not trimming short end as low.limit is outside spectral data range.")
        }
      trim.low <- FALSE
    }
  }

  # check whether we should expand the high end
  high.end <- max(spct, na.rm = TRUE)
  if (trim.high && high.end < high.limit) {
    if (!is.null(fill)) {
      # expand short tail
      high.tail.length <- trunc(high.limit - high.end) + ifelse(use.hinges, 2, 1)
      high.tail.w.length <- seq(from = ifelse(use.hinges, high.end + 1e-12, high.end + 1),
                                to = high.limit,
                                length.out = high.tail.length)
      spct.bottom <- dplyr::data_frame(w.length = high.tail.w.length)
      for (data.col in names.data) {
        spct.bottom[[data.col]] <- fill
      }
      spct <- plyr::rbind.fill(list(spct, spct.bottom))
      spct <- dplyr::as_data_frame(spct)
      setGenericSpct(spct)
      low.end <- max(spct)
    } else {
      # give a warning only if difference is > 0.01 nm
      if (verbose && (high.limit - high.end) > 0.01) {
        warning("Not trimming long end as high.limit is outside spectral data range.")
      }
     trim.high <- FALSE
    }
  }

  # insert hinges
  if (use.hinges) {
    hinges <- NULL
    if (trim.low) {
      hinges <- c(hinges, low.limit - 1e-12, low.limit)
    }
    if (trim.high) {
      hinges <- c(hinges, high.limit - 1e-12, high.limit)
    }
    spct <- insert_spct_hinges(spct, hinges)
  } else {

  }
  if (trim.low && trim.high) {
    within.selector <- with(spct, w.length >= low.limit & w.length < high.limit)
  } else if (trim.low) {
    within.selector <- with(spct, w.length >= low.limit)
  } else if (trim.high) {
    within.selector <- with(spct, w.length < high.limit)
  } else {
    within.selector <- TRUE
  }
  if (is.null(fill)) {
    spct <- spct[within.selector, ]
  } else {
    for (data.col in names.data) {
      spct[!within.selector, data.col] <- fill
    }
  }
  #
  spct <- copy_attributes(x, spct)
  # we now use plyr::rbind.fill which does not remove attributes
  # most of the code below may be redundant!!!
  # class(spct) <- class_spct
  # if (!is.null(comment.spct)) {
  #   comment(spct) <- comment.spct
  # }
  # if (!is.null(time.unit.spct) && !is.na(time.unit.spct)) {
  #   setTimeUnit(spct, time.unit.spct)
  # }
  # if (!is.null(Tfr.type.spct)) {
  #   setTfrType(spct, Tfr.type.spct)
  # }
  # if (!is.null(Rfr.type.spct)) {
  #   setRfrType(spct, Rfr.type.spct)
  # }
  # if (!is.null(when.measured.spct)) {
  #   setWhenMeasured(spct, when.measured.spct)
  # }
  # if (!is.null(where.measured.spct)) {
  #   setWhereMeasured(spct, where.measured.spct)
  # }
  # if (!is.null(what.measured.spct)) {
  #   setWhatMeasured(spct, what.measured.spct)
  # }
  # if (!is.null(instr.settings.spct)) {
  #   setInstrSettings(spct, instr.settings.spct)
  # }
  if (byref && is.name(name)) {
    name <- as.character(name)
    assign(name, spct, parent.frame(), inherits = TRUE)
  }
  check_spct(spct)
  spct
}

#' @rdname trim_spct
#'
#' @param mspct an object of class "generic_mspct"
#'
#' @export
#'
trim_mspct <- function(mspct,
                       range = NULL,
                       low.limit = NULL,
                       high.limit = NULL,
                       use.hinges = TRUE,
                       fill = NULL,
                       byref = FALSE,
                       verbose = getOption("photobiology.verbose", default = TRUE) ) {
  name <- substitute(mspct)

  z <- msmsply(mspct = mspct,
               .fun = trim_spct,
               range = range,
               low.limit = low.limit,
               high.limit = high.limit,
               use.hinges = use.hinges,
               fill = fill,
               byref = FALSE,
               verbose = verbose )

  if (byref & is.name(name)) {
    name <- as.character(name)
    assign(name, z, parent.frame(), inherits = TRUE)
  }
  z
}

#' Trim head and/or tail of a spectrum
#'
#' Triming of head and tail of a spectrum based on wavelength limits,
#' interpolation used by default. Expansion is also possible.
#'
#' @param x an R object
#' @param range a numeric vector of length two, or any other object for which
#'   function range() will return two
#' @param use.hinges logical, if TRUE (the default) wavelengths in nm.
#' @param fill if \code{fill == NULL} then tails are deleted, otherwise tails
#'   are filled with the value of fill.
#' @param ... not used
#'
#' @return an R object of same class as input, usually of a different
#'   length, either shorter or longer.
#'
#' @note By default the \code{w.length} values for the first and last rows
#'   in the returned object are the values supplied as \code{range}.
#'
#' @family trim functions
#' @export
#' @examples
#' trim_wl(sun.spct, range = c(400, 500))
#' trim_wl(sun.spct, range = c(NA, 500))
#' trim_wl(sun.spct, range = c(400, NA))
#'
trim_wl <- function(x, range, use.hinges, fill, ...) UseMethod("trim_wl")

#' @describeIn trim_wl Default for generic function
#'
#' @export
#'
trim_wl.default <- function(x, range, use.hinges, fill, ...) {
  warning("'trim_wl' is not defined for objects of class ", class(x)[1])
  x
}

#' @describeIn trim_wl Trim an object of class "generic_spct" or derived.
#'
#' @export
#'
trim_wl.generic_spct <- function(x,
                                 range = NULL,
                                 use.hinges = TRUE,
                                 fill = NULL, ...) {
  if (is.null(range)) {
    return(x)
  }
  trim_spct(spct = x,
            range = range,
            use.hinges = use.hinges,
            fill = fill,
            byref = FALSE,
            verbose = getOption("photobiology.verbose", default = FALSE) )
}

#' @describeIn trim_wl  Trim an object of class "generic_mspct" or derived.
#'
#' @export
#'
trim_wl.generic_mspct <- function(x,
                                  range = NULL,
                                  use.hinges = TRUE,
                                  fill = NULL, ...) {
  if (is.null(range)) {
    return(x)
  }
  trim_mspct(mspct = x,
             range = range,
             use.hinges = use.hinges,
             fill = fill,
             byref = FALSE,
             verbose = getOption("photobiology.verbose", default = FALSE) )
}

#' @describeIn trim_wl Trim an object of class "waveband".
#' @param trim logical (default is TRUE which trims the wavebands at the
#'   boundary, while FALSE discards wavebands that are partly off-boundary).
#'
#' @note trim_wl when applied to waveband objects always inserts hinges when
#'   trimming.
#'
#' @export
#'
trim_wl.waveband <- function(x,
                             range = NULL,
                             use.hinges = TRUE,
                             fill = NULL,
                             trim = getOption("photobiology.waveband.trim",
                                              default = TRUE),
                             ...) {
  if (is.null(range)) {
    return(x)
  }
  trim_waveband(w.band = x,
            range = range,
            trim = trim,
            use.hinges = use.hinges)
}

#' @describeIn trim_wl Trim a list (of "waveband" objects).
#'
#' @note trim_wl when applied to waveband objects always inserts hinges when
#'   trimming.
#'
#' @export
#'
trim_wl.list <- function(x,
                         range = NULL,
                         use.hinges = TRUE,
                         fill = NULL,
                         trim = getOption("photobiology.waveband.trim",
                                          default = TRUE),
                         ...) {
  stopifnot(all(sapply(x, is.waveband)))
  if (is.null(range)) {
    return(x)
  }
  trim_waveband(w.band = x,
                range = range,
                trim = trim,
                use.hinges = use.hinges)
}

#' Clip head and/or tail of a spectrum
#'
#' Clipping of head and tail of a spectrum based on wavelength limits, no
#' interpolation used.
#'
#' @param x an R object
#' @param range a numeric vector of length two, or any other object for which
#'   function \code{range()} will return range of walengths expressed in
#'   nanometres.
#' @param ... not used
#'
#' @return an R object of same class as input, most frequently of a shorter
#'   length, and never longer.
#'
#' @note The condition tested is \code{wl >= range[1] & wl < (range[2] + 1e-13)}.
#'
#' @family trim functions
#' @export
#' @examples
#' clip_wl(sun.spct, range = c(400, 500))
#' clip_wl(sun.spct, range = c(NA, 500))
#' clip_wl(sun.spct, range = c(400, NA))
#'
clip_wl <- function(x, range, ...) UseMethod("clip_wl")

#' @describeIn clip_wl Default for generic function
#'
#' @export
#'
clip_wl.default <- function(x, range, ...) {
  warning("'clip_wl' is not defined for objects of class ", class(x)[1])
  x
}

#' @describeIn clip_wl Clip an object of class "generic_spct" or derived.
#'
#' @export
#'
clip_wl.generic_spct <- function(x, range = NULL, ...) {
  if (is.null(range)) {
    return(x)
  }
  guard <- 1e-13
  stopifnot(is.any_spct(x))
  stopifnot(!all(is.na(range)))
  if (is.numeric(range) && length(range) == 2) {
    if (is.na(range[1])) {
      x[x[["w.length"]] < range[2] + guard, ]
    } else if (is.na(range[2])) {
      x[x[["w.length"]] >= range[1], ]
    } else {
      x[x[["w.length"]] >= range[1] & x[["w.length"]] < range[2] + guard, ]
    }
  } else {
    range = range(range)
    x[x[["w.length"]] >= range[1] & x[["w.length"]] < range[2] + guard, ]
  }
}

#' @describeIn clip_wl  Clip an object of class "generic_mspct" or derived.
#'
#' @export
#'
clip_wl.generic_mspct <- function(x, range = NULL, ...) {
  msmsply(mspct = x,
          .fun = clip_wl,
          range = range)
}

#' @describeIn clip_wl Clip an object of class "waveband".
#'
#' @export
#'
clip_wl.waveband <- function(x, range = NULL, ...) {
  if (is.null(range)) {
    return(x)
  }
  trim_waveband(w.band = x,
                range = range,
                trim = FALSE)
}

#' @describeIn clip_wl Clip a list (of objects of class "waveband").
#'
#' @export
#'
clip_wl.list <- function(x, range = NULL, ...) {
  stopifnot(all(sapply(x, is.waveband)))
  if (is.null(range)) {
    return(x)
  }
  trim_waveband(w.band = x,
                range = range,
                trim = FALSE)
}
