#' Absorbance
#'
#' Function to calculate the mean, total, or other summary of absorbance for
#' spectral data stored in a \code{filter_spct} or in an \code{object_spct}.
#'
#' @param spct an R object
#' @param w.band waveband or list of waveband objects The waveband(s) determine
#'   the region(s) of the spectrum that are summarized.
#' @param quantity character
#' @param wb.trim logical Flag indicating if wavebands crossing spectral data
#'   boundaries are trimmed or ignored
#' @param use.hinges logical Flag indicating whether to use hinges to reduce
#'   interpolation errors
#' @param ... other arguments (possibly ignored)
#'
#' @note The \code{use.hinges} parameter controls speed optimization. The
#'   defaults should be suitable in most cases. Only the range of wavelengths in
#'   the wavebands is used and all BSWFs are ignored.
#'
#' @examples
#' absorbance(polyester.spct, new_waveband(400,700))
#' absorbance(yellow_gel.spct, new_waveband(400,700))
#'
#' @export
#'
absorbance <- function(spct, w.band, quantity, wb.trim, use.hinges, ...) UseMethod("absorbance")

#' @describeIn absorbance Default for generic function
#'
#' @export
#'
absorbance.default <- function(spct, w.band, quantity, wb.trim, use.hinges, ...) {
  warning("'absorbance' is not defined for objects of class ", class(spct)[1])
  return(NA_real_)
}

#' @describeIn absorbance Specialization for filter spectra
#'
#' @export
#'
absorbance.filter_spct <-
  function(spct, w.band=NULL, quantity="average",
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL), ...) {
    absorbance_spct(spct = spct, w.band = w.band, quantity = quantity, wb.trim = wb.trim,
                    use.hinges = use.hinges)
  }

#' @describeIn absorbance Specialization for object spectra
#'
#' @export
#'
absorbance.object_spct <-
  function(spct, w.band=NULL, quantity="average",
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL), ...) {
    spct <- as.filter_spct(spct)
    absorbance_spct(spct, w.band = w.band, quantity = quantity, wb.trim = wb.trim, use.hinges = use.hinges)
  }

#' Calculate absorbance from spectral absorbance.
#'
#' This function returns the mean absorbance for a given
#' waveband of a absorbance spectrum.
#'
#' @param spct filter_spct
#' @param w.band waveband or list of waveband objects The wavebands determine
#'   the region(s) of the spectrum that are summarized.
#' @param quantity character
#' @param wb.trim logical Flag if wavebands crossing spectral data boundaries
#'   are trimmed or ignored
#' @param use.hinges logical Flag indicating whether to use hinges to reduce
#'   interpolation errors
#' @keywords internal
#'
absorbance_spct <-
  function(spct, w.band, quantity, wb.trim, use.hinges) {
    if (is_normalized(spct) || is_scaled(spct)) {
      warning("The spectral data has been normalized or scaled, making impossible to calculate absorbance")
      return(NA)
    }
    spct <- T2A(spct, action="replace", byref=FALSE)
    spct <- spct[ , c("w.length", "A")]
    # if the waveband is undefined then use all data
    if (is.null(w.band)){
      w.band <- waveband(spct)
    }
    if (is.waveband(w.band)) {
      # if the argument is a single w.band, we enclose it in a list
      # so that the for loop works as expected.This is a bit of a
      # cludge but let's us avoid treating it as a special case
      w.band <- list(w.band)
    }
    w.band <- trim_waveband(w.band=w.band, range=spct, trim=wb.trim)

    # if the w.band includes 'hinges' we insert them
    # choose whether to use hinges or not
    # if the user has specified its value, we leave it alone
    # but if it was not requested, we decide whether to insert
    # hinges or not based of the wavelength resolution of the
    # spectrum. This will produce small errors for high
    # spectral resolution data, and speed up the calculations
    # a lot in such cases
    if (is.null(use.hinges)) {
      use.hinges <- auto_hinges(spct[["w.length"]])
    }

    # we collect all hinges and insert them in one go
    # this may alter a little the returned values
    # but should be faster
    if (use.hinges) {
      all.hinges <- NULL
      for (wb in w.band) {
        if (!is.null(wb$hinges) & length(wb[["hinges"]]>0)) {
          all.hinges <- c(all.hinges, wb[["hinges"]])
        }
      }
      if (!is.null(all.hinges)) {
        spct <- insert_spct_hinges(spct, all.hinges)
      }
    }

    # we prepare labels for output
    wb.name <- names(w.band)
    no_names_flag <- is.null(wb.name)
    if (no_names_flag) {
      wb.name <- character(length(w.band))
    }
    #
    # we iterate through the list of wavebands
    absorbance <- numeric(length(w.band))
    i <- 0
    for (wb in w.band) {
      i <- i + 1
      # we get names from wb if needed
      if (no_names_flag) {
        if (is_effective(wb)) {
          warning("Using only wavelength range from a weighted waveband object.")
          wb.name[i] <- paste("range", as.character(signif(min(wb), 4)),
                              as.character(signif(max(wb), 4)), sep=".")
        } else {
          wb.name[i] <- wb$name
        }
      }
      # we calculate the average transmittance.
      absorbance[i] <- integrate_spct(trim_spct(spct, wb, use.hinges=FALSE))
    }

    if (quantity %in% c("contribution", "contribution.pc")) {
      total <- absorbance_spct(spct, w.band=NULL,
                                  quantity="total", use.hinges=use.hinges)
      absorbance <- absorbance / total
      if (quantity == "contribution.pc") {
        absorbance <- absorbance * 1e2
      }
    } else if (quantity %in% c("relative", "relative.pc")) {
      total <- sum(absorbance)
      absorbance <- absorbance / total
      if (quantity == "relative.pc") {
        absorbance <- absorbance * 1e2
      }
    } else if (quantity %in% c("average", "mean")) {
      absorbance <- absorbance / sapply(w.band, spread)
    }
    if (length(absorbance) == 0) {
      absorbance <- NA
      names(absorbance) <- "out of range"
    }
    names(absorbance) <- paste(names(absorbance), wb.name)
    attr(absorbance, "Tfr.type") <- getTfrType(spct)
    attr(absorbance, "radiation.unit") <- paste("absorbance", quantity)
    return(absorbance)
  }

# filter_mspct methods -----------------------------------------------

#' @describeIn absorbance Calculates absorbance from a \code{filter_mspct}
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
absorbance.filter_mspct <-
  function(spct, w.band=NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = absorbance,
      w.band = w.band,
      quantity = quantity,
      wb.trim = wb.trim,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }

# object_mspct methods -----------------------------------------------

#' @describeIn absorbance Calculates absorbance from a \code{object_mspct}
#'
#' @export
#'
absorbance.object_mspct <-
  function(spct, w.band=NULL,
           quantity="average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges=getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = absorbance,
      w.band = w.band,
      quantity = quantity,
      wb.trim = wb.trim,
      use.hinges = use.hinges,
      col.names = names(w.band),
      idx = idx
    )
  }
