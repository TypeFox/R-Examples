#' Reflectance
#'
#' Function to calculate the mean, total, or other summary of reflectance for
#' spectral data stored in a \code{reflector_spct} or in an \code{object_spct}.
#'
#' @param spct an R object
#' @param w.band waveband or list of waveband objects The waveband(s) determine
#'   the region(s) of the spectrum that are summarized.
#' @param quantity character
#' @param wb.trim logical Flag telling if wavebands crossing spectral data boundaries
#'   are trimmed or ignored
#' @param use.hinges logical Flag indicating whether to use hinges to reduce
#'   interpolation errors
#' @param ... other arguments
#'
#' @note The \code{use.hinges} parameter controls speed optimization. The
#'   defaults should be suitable in most cases. Only the range of wavelengths
#'   in the wavebands is used and all BSWFs are ignored.
#'
#' @return A single numeric value with no change in scale factor
#'
#' @examples
#' reflectance(black_body.spct, waveband(c(400,700)))
#' reflectance(white_body.spct, waveband(c(400,700)))
#'
#' @export
#'
reflectance <- function(spct, w.band, quantity, wb.trim, use.hinges, ...) UseMethod("reflectance")

#' @describeIn reflectance Default for generic function
#'
#' @export
#'
reflectance.default <- function(spct, w.band, quantity, wb.trim, use.hinges, ...) {
  warning("'reflectance' is not defined for objects of class ", class(spct)[1])
  return(NA)
}

#' @describeIn reflectance Specialization for reflector_spct
#'
#' @export
#'
reflectance.reflector_spct <-
  function(spct, w.band = NULL, quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL), ... ) {
    reflectance_spct(spct = spct, w.band = w.band,
                     quantity = quantity,
                     wb.trim = wb.trim,
                     use.hinges = use.hinges)
  }

#' @describeIn reflectance Specialization for object_spct
#'
#' @export
#'
reflectance.object_spct <-
  function(spct, w.band = NULL, quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default =TRUE),
           use.hinges=getOption("photobiology.use.hinges", default=NULL), ... ) {
    reflectance_spct(spct = spct, w.band = w.band,
                     quantity = quantity,
                     wb.trim = wb.trim,
                     use.hinges = use.hinges)
  }

#' Calculate reflectance from spectral reflectance
#'
#' This function returns the mean reflectance for a given waveband and a
#' reflectance spectrum.
#'
#' @param spct an object of class generic_spct"
#' @param w.band list of waveband definitions created with new_waveband()
#' @param quantity character string
#' @param wb.trim logical if TRUE wavebands crossing spectral data boundaries
#'   are trimmed, if FALSE, they are discarded
#' @param use.hinges logical indicating whether to use hinges to reduce
#'   interpolation errors
#'
#' @return A single numeric value expressed as a fraction of one
#' @keywords internal
#'
reflectance_spct <-
  function(spct, w.band, quantity, wb.trim, use.hinges){
    Rfr.type <- getRfrType(spct)
    if (is.object_spct(spct)) {
      spct <- as.reflector_spct(spct)
    }
    spct <- spct[ , c("w.length", "Rfr")]
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
        if (!is.null(wb$hinges) && length(wb$hinges) > 0) {
          all.hinges <- c(all.hinges, wb$hinges)
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
   # we iterate through the list of wavebands
    reflectance <- numeric(length(w.band))
    i <- 0
    for (wb in w.band) {
      i <- i + 1
      # we get names from wb if needed
      if (no_names_flag) {
        if (is_effective(wb)) {
          warning("Using only wavelength range from a weighted waveband object.")
          wb.name[i] <- paste("range",
                              as.character(signif(min(wb), 4)),
                              as.character(signif(max(wb), 4)), sep=".")
        } else {
          wb.name[i] <- wb$name
        }
      }
      # we calculate the average reflectance.
      reflectance[i] <- integrate_spct(trim_spct(spct, wb, use.hinges = FALSE))
    }

   if (quantity %in% c("contribution", "contribution.pc")) {
     total <- reflectance_spct(spct, w.band = NULL, wb.trim = wb.trim,
                                quantity = "total", use.hinges = use.hinges)
     reflectance <- reflectance / total
     if (quantity == "contribution.pc") {
       reflectance <- reflectance * 1e2
     }
   } else if (quantity %in% c("relative", "relative.pc")) {
     total <- sum(reflectance)
     reflectance <- reflectance / total
     if (quantity == "relative.pc") {
       reflectance <- reflectance * 1e2
     }
   } else if (quantity %in% c("average", "mean")) {
     reflectance <- reflectance / sapply(w.band, spread)
   } else if (quantity == "total") {
   } else if (quantity != "total") {
     warning("'quantity '", quantity, "' is invalid, returning 'total' instead")
     quantity <- "total"
   }
   if (length(reflectance) == 0) {
     reflectance <- NA
     names(reflectance) <- "off range"
   }
   names(reflectance) <- paste(names(reflectance), wb.name)
   attr(reflectance, "Rfr.type") <- getRfrType(spct)
   attr(reflectance, "radiation.unit") <- paste("reflectance", quantity)
   return(reflectance)
  }

# reflector_mspct methods -----------------------------------------------

#' @describeIn reflectance Calculates reflectance from a \code{reflector_mspct}
#'
#' @param idx logical whether to add a column with the names of the elements of
#'   spct
#'
#' @export
#'
reflectance.reflector_mspct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges = getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = reflectance,
      w.band = w.band,
      quantity = quantity,
      wb.trim = wb.trim,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }

# object_mspct methods -----------------------------------------------

#' @describeIn reflectance Calculates reflectance from a \code{object_mspct}
#'
#' @export
#'
reflectance.object_mspct <-
  function(spct, w.band = NULL,
           quantity = "average",
           wb.trim = getOption("photobiology.waveband.trim", default = TRUE),
           use.hinges= getOption("photobiology.use.hinges", default = NULL),
           ..., idx = !is.null(names(spct))) {
    msdply(
      mspct = spct,
      .fun = reflectance,
      w.band = w.band,
      quantity = quantity,
      wb.trim = wb.trim,
      use.hinges = use.hinges,
      idx = idx,
      col.names = names(w.band)
    )
  }
