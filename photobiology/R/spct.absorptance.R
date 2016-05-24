#' Absorptance
#'
#' Function to calculate the mean, total, or other summary of absorptance for
#' spectral data stored in a \code{filter_spct} or in an \code{object_spct}.
#' Absorptance is a different quantity than absorbance.
#'
#' @param spct an R object
#' @param w.band waveband or list of waveband objects The waveband(s) determine
#'   the region(s) of the spectrum that are summarized.
#' @param quantity character
#' @param wb.trim logical Flag if wavebands crossing spectral data boundaries
#'   are trimmed or ignored
#' @param use.hinges logical Flag indicating whether to use hinges to reduce
#'   interpolation errors
#' @param ... other arguments (possibly ignored)
#'
#' @note The \code{use.hinges} parameter controls speed optimization. The
#'   defaults should be suitable in most cases. Only the range of wavelengths
#'   in the wavebands is used and all BSWFs are ignored.
#'
#' @return A single numeric value with no change in scale factor, except in the
#' case of percentages (absorptance is the fraction absorbed)
#'
#' @examples
#' absorptance(black_body.spct, new_waveband(400,500))
#' absorptance(white_body.spct, new_waveband(300,400))
#'
#' @export
#'
absorptance <- function(spct, w.band, quantity, wb.trim, use.hinges, ...) UseMethod("absorptance")

#' @describeIn absorptance Default for generic function
#'
#' @export
#'
absorptance.default <- function(spct, w.band, quantity, wb.trim, use.hinges, ...) {
  warning("'absorptance' is not defined for objects of class ", class(spct)[1])
  return(NA_real_)
}

#' @describeIn absorptance Specialization for filter spectra
#'
#' @export
#'
absorptance.filter_spct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL), ... ) {
    if (getTfrType(spct) != "internal") {
      warning("Internal absorptance cannot be calculed from total transmittance alone")
      return(NA)
    } else {
      absorptance_spct(spct = spct, w.band = w.band, quantity = quantity,
                       wb.trim = wb.trim, use.hinges = use.hinges)
    }
  }

#' @describeIn absorptance Specialization for object spectra
#'
#' @export
#'
absorptance.object_spct <-
  function(spct, w.band=NULL,
           quantity="average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges=getOption("photobiology.use.hinges", default = NULL), ...)  {
    absorptance_spct(spct = spct, w.band = w.band, quantity = quantity,
                     wb.trim = wb.trim, use.hinges = use.hinges)
  }

#' Calculate absorptance from spectral absorptance.
#'
#' This function returns the summary absorptance for a given waveband of a
#' \code{object_spct} object
#'
#' @param spct object_spct
#' @param w.band waveband or list of waveband objects The wavebands determine
#'   the region(s) of the spectrum that are summarized.
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.hinges logical indicating whether to use hinges to reduce
#'   interpolation errors
#'
#' @keywords internal
#'
absorptance_spct <-
  function(spct, w.band = NULL, quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL) ) {
    if (is_normalized(spct) || is_scaled(spct)) {
      warning("The spectral data has been normalized or scaled, making impossible to calculate absorptance")
      return(NA)
    }

    # we calculate absorptance
    Tfr.type <- getTfrType(spct)
    Rfr.type <- getRfrType(spct)
    if (is.filter_spct(spct) && Tfr.type == "internal") {
      Afr.type <- Tfr.type
      Rfr.type <- "unknown" # otherwise NA would require special handling
      A2T(spct, action = "add", byref = TRUE)
      temp.spct <- dplyr::data_frame(w.length = spct[["w.length"]],
                                     Afr = 1 - spct[["Tfr"]])
    } else if (Tfr.type == "total" && Rfr.type == "total") {
      Afr.type <- "total"
      temp.spct <- dplyr::data_frame(w.length = spct[["w.length"]],
                               Afr = 1 - spct[["Tfr"]] - spct[["Rfr"]])
     } else if (Tfr.type == "internal" && Rfr.type == "total") {
      Afr.type <- "total"
      temp.spct <- dplyr::data_frame(w.length = spct[["w.length"]],
                                   Afr = (1 - spct[["Tfr"]]) * (1 - spct[["Rfr"]]))
    } else if (Tfr.type == "unknown" || Rfr.type == "unknown") {
      warning("'unknown' Tfr.type or Rfr.type, skipping absorptance calculation")
      absorptance <- NA
      attr(absorptance, "Afr.type") <- "unknown"
      attr(absorptance, "radiation.unit") <- paste("absorptance", quantity)
      return(absorptance)
    } else if (Rfr.type == "specular") {
      warning("'specular' Rfr.type, skipping absorptance calculation")
      absorptance <- NA
      attr(absorptance, "Afr.type") <- "unknown"
      attr(absorptance, "radiation.unit") <- paste("absorptance", quantity)
      return(absorptance)
    } else {
      stop("Failed assertion with Tfr.type: ", Tfr.type, "and Rfr.type: ", Rfr.type)
    }
    temp.spct <- setGenericSpct(temp.spct)
    # if the waveband is undefined then use all data
    if (is.null(w.band)) {
      w.band <- waveband(spct)
    }
    if (is.waveband(w.band)) {
      # if the argument is a single w.band, we enclose it in a list
      # so that the for loop works as expected.This is a bit of a
      # cludge but let's us avoid treating it as a special case
      w.band <- list(w.band)
    }
    w.band <- trim_waveband(w.band = w.band, range = spct, trim = wb.trim)

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
        if (!is.null(wb$hinges) & length(wb$hinges) > 0) {
          all.hinges <- c(all.hinges, wb$hinges)
        }
      }
      if (!is.null(all.hinges)) {
        temp.spct <- insert_spct_hinges(temp.spct, all.hinges)
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
    absorptance <- numeric(length(w.band))
    i <- 0
    for (wb in w.band) {
      i <- i + 1
      # we get names from wb if needed
      if (no_names_flag) {
        if (is_effective(wb)) {
          warning("Using only wavelength range from a weighted waveband object.")
          wb.name[i] <- paste("range", as.character(signif(min(wb), 4)), as.character(signif(max(wb), 4)), sep = ".")
        } else {
          wb.name[i] <- wb$name
        }
      }
      absorptance[i] <-
        integrate_spct(trim_spct(temp.spct, wb,
                                 use.hinges = FALSE))
    }

    if (quantity %in% c("contribution", "contribution.pc")) {
      total <- absorptance_spct(spct, w.band = NULL,
                                quantity = "total",
                                use.hinges = use.hinges)
      absorptance <- absorptance / total
      if (quantity == "contribution.pc") {
        absorptance <- absorptance * 1e2
      }
    } else if (quantity %in% c("relative", "relative.pc")) {
      total <- sum(absorptance)
      absorptance <- absorptance / total
      if (quantity == "relative.pc") {
        absorptance <- absorptance * 1e2
      }
    } else if (quantity %in% c("average", "mean")) {
      absorptance <- absorptance / sapply(w.band, spread)
    }
    if (length(absorptance) == 0) {
      absorptance <- NA
      names(absorptance) <- "out of range"
    }
    names(absorptance) <- paste(names(absorptance), wb.name)
    attr(absorptance, "Afr.type") <- Afr.type
    attr(absorptance, "radiation.unit") <- paste("absorptance", quantity)
    return(absorptance)
  }

#' @describeIn absorptance Calculates absorptance from a \code{filter_mspct}
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
absorptance.filter_mspct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct)) ) {
    msdply(
      mspct = spct,
      .fun = absorptance,
      w.band = w.band,
      quantity = quantity,
      wb.trim = wb.trim,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }

#' @describeIn absorptance Calculates absorptance from a \code{object_mspct}
#'
#' @export
#'
absorptance.object_mspct <-
  function(spct, w.band=NULL,
           quantity="average",
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL),
           ..., idx = !is.null(names(spct)) ) {
    msdply(
      mspct = spct,
      .fun = absorptance,
      w.band = w.band,
      quantity = quantity,
      wb.trim = wb.trim,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }
