# response methods --------------------------------------------------------


#' Integrated response
#'
#' Calculate average photon- or energy-based photo-response.
#'
#' @param spct an R object of class "generic_spct"
#' @param w.band waveband or list of waveband objects The waveband(s) determine
#'   the region(s) of the spectrum that are summarized
#' @param unit.out character Allowed values "energy", and "photon", or its alias
#'   "quantum"
#' @param quantity character Allowed values ""
#' @param time.unit character or lubridate::duration
#' @param wb.trim logical Flag telling if wavebands crossing spectral data boundaries
#'   are trimmed or ignored
#' @param use.hinges logical indicating whether to use hinges to reduce
#'   interpolation errors
#' @param ... other arguments
#'
#' @note The parameter \code{use.hinges} controls speed optimization. The
#'   defaults should be suitable in mosts cases. Only the range of wavelengths
#'   in the wavebands is used and all BSWFs are ignored.
#'
#' @return A single numeric value expressed either as a fraction of one or a
#'   percentage, or a vector of the same length as the list of wave.bands. The
#'   quantity returned depends on the value of \code{quantity}. Whether it is
#'   expressed in energy-based or photon-based units depends on \code{unit.out}.
#'
#' @export
#' @family response functions
#'
response <- function(spct, w.band, unit.out, quantity, time.unit, wb.trim, use.hinges, ...) UseMethod("response")

#' @describeIn response Default for generic function
#'
#' @export
#'
response.default <- function(spct, w.band, unit.out, quantity, time.unit, wb.trim, use.hinges, ...) {
  warning("'response' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn response Method for response spectra.
#'
#' @export
#'
response.response_spct <-
  function(spct, w.band = NULL,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           quantity = "total",
           time.unit = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL), ... ) {
    resp_spct(spct = spct, w.band = w.band, unit.out = unit.out,
              quantity = quantity, time.unit = time.unit, wb.trim = wb.trim,
              use.hinges = use.hinges )
  }

#' Calculate response from spectral response
#'
#' This function returns the mean response for a given waveband and a response
#' spectrum.
#'
#' @param spct an object of class response_spct"
#' @param w.band a waveband object or a list of waveband objects
#' @param unit.out character with allowed values "energy", and "photon", or its
#'   alias "quantum"
#' @param quantity character with allowed values "total", "average" ("mean"),
#'   "contibution", "contribution.pc", "relative", "relative.pc"
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.hinges logical indicating whether to use hinges to reduce
#'   interpolation errors
#' @param ... other arguments
#'
#' @return a single numeric value expressed either as a fraction of one or a
#'   percentage, or a vector of the same length as the list of wave.bands.
#' @keywords internal
#'
#' @note The parameter \code{use.hinges} controls speed optimization. The
#'   defaults should be suitable in mosts cases. Only the range of wavelengths
#'   in the wavebands is used and all BSWFs are ignored.
#'
#' @keywords internal
#'
resp_spct <-
  function(spct, w.band, unit.out, quantity, time.unit, wb.trim, use.hinges, ...) {
    if (is_normalized(spct) || is_scaled(spct)) {
      warning("The spectral data has been normalized or scaled, ",
              "making impossible to calculate integrated response")
      return(NA)
    }
    # makes "quantum" synonym for "photon" without changes to other code
    if (unit.out == "quantum") {
      unit.out <- "photon"
    }

    data.time.unit <- getTimeUnit(spct, force.duration = lubridate::is.duration(time.unit))

    if (!is.null(time.unit) && time.unit != data.time.unit) {
      if (!lubridate::is.duration(time.unit) && !is.character(time.unit)) {
        message("converting 'time.unit' ", time.unit, " into a lubridate::duration")
        time.unit <- lubridate::as.duration(time.unit)
      }
      spct <- convertTimeUnit(spct, time.unit = time.unit, byref = FALSE)
    } else {
      time.unit <- data.time.unit
    }

    if (unit.out=="photon") {
      spct <- e2q(spct)
      spct <- spct[ , c("w.length", "s.q.response")]
    } else if (unit.out=="energy") {
      spct <- q2e(spct)
      spct <- spct[ , c("w.length", "s.e.response")]
    } else {
      stop("Invalid 'unit.out'")
    }

    # if the waveband is undefined then use all data
    if (is.null(w.band)){
      w.band <- waveband(spct)
    }
    if (is.waveband(w.band)) {
      # if the argument is a single w.band, we enclose it in a list
      # so that the for loop works as expected. This is a bit of a
      # cludge but it let's us avoid treating it as a special case
      w.band <- list(w.band)
    }
    w.band <- trim_waveband(w.band=w.band, range=spct, trim=wb.trim)

    # if the w.band includes 'hinges' we insert them,
    # but if not, we decide whether to insert hinges or not
    # hinges or not based of the wavelength resolution of the
    # spectrum. This can produce small errors for high
    # spectral resolution data, but speed up the calculations.
    if (is.null(use.hinges)) {
      use.hinges <- auto_hinges(spct[["w.length"]])
    }

    # we collect all hinges and insert them in one go
    # this may alter very slightly the returned values
    # but improves calculation speed
    if (use.hinges) {
      all.hinges <- NULL
      for (wb in w.band) {
        if (!is.null(wb$hinges) && length(wb$hinges)>0) {
          all.hinges <- c(all.hinges, wb$hinges)
        }
      }
      if (!is.null(all.hinges)) {
        spct <- insert_spct_hinges(spct, all.hinges)
      }
    }

    # we prepare labels for output
    wb_name <- names(w.band)
    no_names_flag <- is.null(wb_name)
    if (no_names_flag) {
      wb_name <- character(length(w.band))
    }

    # we iterate through the list of wavebands
    response <- double(length(w.band))
    i <- 0
    for (wb in w.band) {
      i <- i + 1
      # we get names from wb if needed
      if (no_names_flag) {
        if (is_effective(wb)) {
          warning("Using only wavelength range from a weighted waveband object.")
          wb_name[i] <- paste("range", as.character(signif(min(wb), 4)),
                              as.character(signif(max(wb), 4)), sep=".")
        } else {
          wb_name[i] <- wb$name
        }
      }
      # we calculate the integrated response.
      response[i] <- integrate_spct(trim_spct(spct, wb, use.hinges = FALSE))
    }
    if (quantity %in% c("contribution", "contribution.pc")) {
      if (any(sapply(w.band, is_effective))) {
        warning("'quantity '", quantity,
                "' not supported when using BSWFs, returning 'total' instead")
        quantity <- "total"
      } else {
        total <- resp_spct(spct, w.band = NULL, unit.out = unit.out,
                           quantity = "total", time.unit = time.unit,
                           wb.trim = FALSE, use.hinges = use.hinges)
        response <- response / total
        if (quantity == "contribution.pc") {
          response <- response * 1e2
        }
      }
    } else if (quantity %in% c("relative", "relative.pc")) {
      if (any(sapply(w.band, is_effective))) {
        warning("'quantity '", quantity,
                "' not supported when using BSWFs, returning 'total' instead")
        quantity <- "total"
      } else {
        total <- sum(response)
        response <- response / total
        if (quantity == "relative.pc") {
          response <- response * 1e2
        }
      }
    } else if (quantity %in% c("average", "mean")) {
      response <- response / sapply(w.band, spread)
    } else if (quantity != "total") {
      warning("'quantity '", quantity, "' is invalid, returning 'total' instead")
      quantity <- "total"
    }

    if (length(response) == 0) {
      response <- NA
      names(response) <- "out of range"
    }
    names(response) <- paste(names(response), wb_name)
    attr(response, "time.unit") <- getTimeUnit(spct)
    attr(response, "radiation.unit") <- paste(unit.out, "response", quantity)
    return(response)
  }

# e_response methods --------------------------------------------------------

#' Energy-based photo-response
#'
#' This function returns the mean, total, or contribution of response for each
#' waveband and a response spectrum.
#'
#' @param spct an R object
#' @param w.band a waveband object or a list of waveband objects
#' @param quantity character
#' @param time.unit character or lubridate::duration
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.hinges logical indicating whether to use hinges to reduce
#'   interpolation errors
#' @param ... other arguments
#'
#' @return A single numeric value expressed either as a fraction of one or a
#'   percentage, or a vector of the same length as the list of wave.bands. The
#'   quantity returned, although always on energy-based units, depends on the
#'   value of \code{quantity}.
#'
#' @export
#' @examples
#' e_response(ccd.spct, new_waveband(200,300))
#' e_response(photodiode.spct)
#'
#' @note The parameter \code{use.hinges} controls speed optimization. The
#'   defaults should be suitable in mosts cases. Only the range of wavelengths
#'   in the wavebands is used and all BSWFs are ignored.
#'
#' @family response functions
#'
e_response <- function(spct, w.band, quantity, time.unit, wb.trim, use.hinges, ...) UseMethod("e_response")

#' @describeIn e_response Default method for generic function
#'
#' @export
#'
e_response.default <- function(spct, w.band, quantity, time.unit, wb.trim, use.hinges, ...) {
  warning("'e_response' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn e_response Method for response spectra.
#'
#' @export
#'
e_response.response_spct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL), ...) {
    resp_spct(spct=spct, w.band=w.band, unit.out="energy",
              quantity=quantity, time.unit=time.unit, wb.trim=wb.trim,
              use.hinges=use.hinges )
  }

# q_response methods --------------------------------------------------------

##' Photon-based photo-response
#'
#' This function returns the mean response for a given
#' waveband and a response spectrum.
#'
#' @param spct an R object
#' @param w.band a waveband object or a list of waveband objects
#' @param quantity character
#' @param time.unit character or lubridate::duration
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.hinges logical indicating whether to use hinges to reduce
#'   interpolation errors
#' @param ... other arguments
#'
#' @return A single numeric value expressed either as a fraction of one or a
#'   percentage, or a vector of the same length as the list of wave.bands. The
#'   quantity returned, although always on photon-based units, depends on the
#'   value of \code{quantity}.
#'
#' @export
#' @examples
#' q_response(ccd.spct, new_waveband(200,300))
#' q_response(photodiode.spct)
#'
#' @note The parameter \code{use.hinges} controls speed optimization. The
#'   defaults should be suitable in mosts cases. Only the range of wavelengths
#'   in the wavebands is used and all BSWFs are ignored.
#'
#' @family response functions
#'
q_response <- function(spct,
                       w.band,
                       quantity,
                       time.unit,
                       wb.trim,
                       use.hinges,
                       ...) UseMethod("q_response")

#' @describeIn q_response Default method for generic function
#'
#' @export
#'
q_response.default <- function(spct, w.band, quantity, time.unit, wb.trim, use.hinges, ...) {
  warning("'q_response' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn q_response Method for response spectra.
#'
#' @export
#'
q_response.response_spct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL), ... ) {
    resp_spct(spct = spct, w.band = w.band, unit.out = "photon",
              quantity = quantity, time.unit = time.unit, wb.trim = wb.trim,
              use.hinges = use.hinges )
  }

# response_mspct methods -----------------------------------------------

#' @describeIn response Calculates response from a \code{response_mspct}
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
response.response_mspct <-
  function(spct, w.band = NULL,
           unit.out = getOption("photobiology.radiation.unit", default = "energy"),
           quantity = "total",
           time.unit = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = response,
      w.band = w.band,
      unit.out = unit.out,
      quantity = quantity,
      time.unit = time.unit,
      wb.trim = wb.trim,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }

#' @describeIn q_response Calculates photon (quantum) response from a
#'   \code{response_mspct}
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
q_response.response_mspct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = q_response,
      w.band = w.band,
      quantity = quantity,
      time.unit = time.unit,
      wb.trim = wb.trim,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }

#' @describeIn e_response Calculates energy response from a
#'   \code{response_mspct}
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
e_response.response_mspct <-
  function(spct, w.band = NULL,
           quantity = "total",
           time.unit = NULL,
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = e_response,
      w.band = w.band,
      quantity = quantity,
      time.unit = time.unit,
      wb.trim = wb.trim,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }
