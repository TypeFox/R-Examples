#' Blank correction
#'
#' @description This function is used to remove blank from eems which can help to reduce the
#' effect of scatter bands.
#'
#' @template template_eem
#' @template template_blank
#' @template template_details_automatic_blank
#'
#' @references Murphy, K. R., Stedmon, C. a., Graeber, D., & Bro, R. (2013).
#'   Fluorescence spectroscopy and multi-way techniques. PARAFAC. Analytical
#'   Methods, 5(23), 6557. https://doi.org/10.1039/c3ay41160e
#'
#'   \url{http://xlink.rsc.org/?DOI=c3ay41160e}
#'
#' @importFrom rlist list.apply list.group list.ungroup
#' @export
#' @examples
#'
#' ## Example 1
#'
#' # Open the fluorescence eem
#' file <- system.file("extdata/cary/scans_day_1", "sample1.csv", package = "eemR")
#' eem <- eem_read(file)
#'
#' plot(eem)
#'
#' # Open the blank eem
#' file <- system.file("extdata/cary/scans_day_1", "nano.csv", package = "eemR")
#' blank <- eem_read(file)
#'
#' plot(blank)
#'
#' # Remove the blank
#' eem <- eem_remove_blank(eem, blank)
#'
#' plot(eem)
#'
#' ## Example 2
#'
#' # Open the fluorescence eem
#' folder <- system.file("extdata/cary/scans_day_1", package = "eemR")
#' eems <- eem_read(folder)
#'
#' plot(eems, which = 3)
#'
#' # Open the blank eem
#' file <- system.file("extdata/cary/scans_day_1", "nano.csv", package = "eemR")
#' blank <- eem_read(file)
#'
#' plot(blank)
#'
#' # Remove the blank
#' eem <- eem_remove_blank(eems, blank)
#'
#' plot(eems, which = 3)
#'
#' # Automatic correction
#' folder <- system.file("extdata/cary/", package = "eemR")
#'
#' # Look at the folder structure
#' list.files(folder, "*.csv", recursive = TRUE)
#'
#' eems <- eem_read(folder, recursive = TRUE)
#' res <- eem_remove_blank(eems)

eem_remove_blank <- function(eem, blank = NA) {

  stopifnot(.is_eemlist(eem) | .is_eem(eem),
            .is_eemlist(blank) | is.na(blank))

  if(is.na(blank)){

    t <- list.group(eem, ~location)
    t <- lapply(t, function(x){class(x) <- "eemlist"; return(x)})

    res <- list.apply(t, .eem_remove_blank)
    res <- list.ungroup(res)
    class(res) <- "eemlist"
    return(res)

  } else {
    .eem_remove_blank(eem, blank)
  }
}

.eem_remove_blank <- function(eem, blank = NA) {

  stopifnot(.is_eemlist(eem) | .is_eem(eem),
            .is_eemlist(blank) | is.na(blank))

  ## It is a list of eems, then call lapply
  if(.is_eemlist(eem)){

    blank_names <- c("nano", "miliq", "milliq", "mq", "blank")

    # if blank is NA then try to split the eemlist into blank and eems
    if(is.na(blank)){
      blank <- eem_extract(eem, blank_names, remove = FALSE, ignore_case = TRUE,
                           verbose = FALSE)
      eem <- eem_extract(eem, blank_names, remove = TRUE, ignore_case = TRUE,
                         verbose = FALSE)

      if(length(blank) != 1 | length(eem) < 1){
        stop("Cannot find blank for automatic correction.", call. = FALSE)
      }
    }

    res <- lapply(eem,
                  eem_remove_blank,
                  blank = blank)

    class(res) <- class(eem)
    return(res)
  }

  #---------------------------------------------------------------------
  # Do the blank subtraction.
  #---------------------------------------------------------------------
  blank <- unlist(blank, recursive = FALSE)
  x <- eem$x - blank$x

  ## Construct an eem object.
  res <- eem(file = eem$sample,
             x = x,
             ex = eem$ex,
             em = eem$em)

  attributes(res) <- attributes(eem)
  attr(res, "is_blank_corrected") <- TRUE

  return(res)
}

#' Remove Raman and Rayleigh scattering
#'
#' @template template_eem
#'
#' @param type A string, either "raman" or "rayleigh".
#' @param order A integer number, either 1 (first order) or 2 (second order).
#' @param width Slit width in nm for the cut. Default is 10 nm.
#'
#' @references
#'
#' Lakowicz, J. R. (2006). Principles of Fluorescence Spectroscopy.
#' Boston, MA: Springer US.#'
#'
#' \url{https://doi.org/10.1007/978-0-387-46312-4}
#'
#' Murphy, K. R., Stedmon, C. a., Graeber, D., & Bro, R. (2013).
#' Fluorescence spectroscopy and multi-way techniques. PARAFAC. Analytical
#' Methods, 5(23), 6557. https://doi.org/10.1039/c3ay41160e#'
#'
#'  \url{http://xlink.rsc.org/?DOI=c3ay41160e}
#'
#' @export
#' @examples
#' # Open the fluorescence eem
#' file <- system.file("extdata/cary/scans_day_1", "sample1.csv", package = "eemR")
#' eem <- eem_read(file)
#'
#' plot(eem)
#'
#' # Remove the scattering
#' eem <- eem_remove_scattering(eem = eem, type = "raman", order = 1, width = 10)
#'
#' plot(eem)

eem_remove_scattering <- function(eem, type, order = 1, width = 10){

  stopifnot(.is_eemlist(eem) | .is_eem(eem),
            all(type %in% c("raman", "rayleigh")),
            is.numeric(order),
            is.numeric(width),
            length(order) == 1,
            length(type) == 1,
            length(width) == 1,
            is_between(order, 1, 2),
            is_between(width, 0, 100))

  ## It is a list of eems, then call lapply
  if(.is_eemlist(eem)){

    res <- lapply(eem,
                  eem_remove_scattering,
                  type = type,
                  order = order,
                  width = width)

    class(res) <- class(eem)
    return(res)
  }

  #---------------------------------------------------------------------
  # Remove the scattering.
  #---------------------------------------------------------------------

  x <- eem$x
  em <- eem$em
  ex <- eem$ex

  if(type == "raman"){
    ex <- .find_raman_peaks(eem$ex)
  }

  ind1 <- mapply(function(x)em <= x, order * ex - width)
  ind2 <- mapply(function(x)em <= x, order * ex + width)

  ind3 <- ifelse(ind1 + ind2 == 1, NA, 1)

  x <- x * ind3

  ## Construct an eem object.
  res <- eem
  res$x <- x

  attributes(res) <- attributes(eem)
  attr(res, "is_scatter_corrected") <- TRUE

  class(res) <- class(eem)

  return(res)
}

.find_raman_peaks <- function(ex){

  # For water, the Raman peak appears at a wavenumber 3600 cm lower than the
  # incident wavenumber. For excitation at 280 nm, the Raman peak from water
  # occurs at 311 nm. Source : Principles of Fluorescence Spectroscopy (2006) -
  # Third Edition.pdf

  ## Convert wavenumber from nm to cm
  ex_wave_number = 1 / ex

  ## For water. 3600 nm = 0.00036 cm
  raman_peaks = ex_wave_number - 0.00036 # I think Stedmon use 3400 TODO

  ## Bring back to nm
  raman_peaks = 1 / raman_peaks

  #raman_peaks <- -(ex / (0.00036 * ex - 1))

  return(raman_peaks)
}

#' Fluorescence Intensity Calibration Using the Raman Scatter Peak of Water
#'
#' @template template_eem
#' @template template_blank
#' @template template_details_automatic_blank
#'
#' @description Normalize fluorescence intensities to the standard scale of
#'   Raman Units (R.U).
#'
#' @details The normalization procedure consists in dividing all fluorescence
#'   intensities by the area (integral) of the Raman peak. The peak is located
#'   at excitation of 350 nm. (ex = 370) betwen 371 nm. and 428 nm in emission
#'   (371 <= em <= 428).
#'
#' @references
#'
#' Lawaetz, A. J., & Stedmon, C. A. (2009). Fluorescence Intensity Calibration
#' Using the Raman Scatter Peak of Water. Applied Spectroscopy, 63(8), 936-940.
#'
#' \url{https://doi.org/10.1366/000370209788964548}
#'
#' Murphy, K. R., Stedmon, C. a., Graeber, D., & Bro, R. (2013). Fluorescence
#' spectroscopy and multi-way techniques. PARAFAC. Analytical Methods, 5(23),
#' 6557.
#'
#' \url{http://xlink.rsc.org/?DOI=c3ay41160e}
#'
#' @return An object of class \code{eem} containing: \itemize{ \item sample The
#'   file name of the eem. \item x A matrix with fluorescence values. \item em
#'   Emission vector of wavelengths. \item ex Excitation vector of wavelengths.
#'   }
#'
#' @export
#' @examples
#' # Open the fluorescence eem
#' file <- system.file("extdata/cary/scans_day_1", "sample1.csv", package = "eemR")
#' eem <- eem_read(file)
#'
#' plot(eem)
#'
#' # Open the blank eem
#' file <- system.file("extdata/cary/scans_day_1", "nano.csv", package = "eemR")
#' blank <- eem_read(file)
#'
#' # Do the normalisation
#' eem <- eem_raman_normalisation(eem, blank)
#'
#' plot(eem)

eem_raman_normalisation <- function(eem, blank = NA) {

  stopifnot(.is_eemlist(eem) | .is_eem(eem),
            .is_eemlist(blank) | is.na(blank))

  if(is.na(blank)){

    t <- list.group(eem, ~location)
    t <- lapply(t, function(x){class(x) <- "eemlist"; return(x)})

    res <- list.apply(t, .eem_raman_normalisation)
    res <- list.ungroup(res)
    class(res) <- "eemlist"
    return(res)

  } else {
    .eem_raman_normalisation(eem, blank)
  }

}

.eem_raman_normalisation <- function(eem, blank = NA){

  stopifnot(.is_eemlist(eem) | .is_eem(eem),
            .is_eemlist(blank) | is.na(blank))

  ## It is a list of eems, then call lapply
  if(.is_eemlist(eem)){

    blank_names <- c("nano", "miliq", "milliq", "mq", "blank")

    # if blank is NA then try to split the eemlist into blank and eems
    if(is.na(blank)){

      blank <- eem_extract(eem, blank_names, remove = FALSE, ignore_case = TRUE,
                           verbose = FALSE)
      eem <- eem_extract(eem, blank_names, remove = TRUE, ignore_case = TRUE,
                         verbose = FALSE)

      if(length(blank) != 1 | length(eem) < 1){
        stop("Cannot find blank for automatic correction.", call. = FALSE)
      }
    }

    res <- lapply(eem,
                  .eem_raman_normalisation,
                  blank = blank)

    class(res) <- class(eem)
    return(res)
  }

  #---------------------------------------------------------------------
  # Do the normalisation.
  #---------------------------------------------------------------------
  blank <- unlist(blank, recursive = FALSE)
  index_ex <- which(blank$ex == 350)
  index_em <- which(blank$em >= 371 & blank$em <= 428)

  x <- blank$em[index_em]
  y <- blank$x[index_em, index_ex]

  if(any(is.na(x)) | any(is.na(y))){
    stop("NA values found in the blank sample. Maybe you removed scattering too soon?", call. = FALSE)
  }

  area <- sum(diff(x) * (y[-length(y)] + y[-1]) / 2)

  cat("Raman area:", area, "\n")

  x <- eem$x / area

  ## Construct an eem object.
  res <- eem(file = eem$sample,
             x = x,
             ex = eem$ex,
             em = eem$em)

  attributes(res) <- attributes(eem)
  attr(res, "is_raman_normalized") <- TRUE

  class(res) <- class(eem)

  return(res)
}

#' Inner-filter effect correction
#'
#' @template template_eem
#'
#' @param pathlength A numeric value indicating the pathlength (in cm) of the
#'   cuvette used for fluorescence measurement. Default is 1 (1cm).
#'
#' @param absorbance A data frame with:
#'
#'   \describe{ \item{wavelength}{A numeric vector containing wavelenghts.}
#'   \item{...}{One or more numeric vectors containing absorbance spectra.}}
#'
#' @section Names matching:
#'
#'   The names of \code{absorbance} variables are expected to match those of the
#'   eems. If the appropriate absorbance spectrum is not found, an uncorrected
#'   eem will be returned and a warning message will be printed.
#'
#' @section Sample dilution:
#'
#'   Kothawala et al. 2013 have shown that a 2-fold dilution was requiered for
#'   sample presenting total absorbance > 1.5. Accordingly, a message will warn
#'   the user if total absorbance is greater than this threshold.
#'
#' @references Parker, C. a., & Barnes, W. J. (1957). Some experiments with
#'   spectrofluorimeters and filter fluorimeters. The Analyst, 82(978), 606.
#'   \url{https://doi.org/10.1039/an9578200606}
#'
#'   Kothawala, D. N., Murphy, K. R., Stedmon, C. A., Weyhenmeyer, G. A., &
#'   Tranvik, L. J. (2013). Inner filter correction of dissolved organic matter
#'   fluorescence. Limnology and Oceanography: Methods, 11(12), 616-630.
#'   \url{https://doi.org/10.4319/lom.2013.11.616}
#'
#' @return An object of class \code{eem} containing: \itemize{ \item sample The
#'   file name of the eem. \item x A matrix with fluorescence values. \item em
#'   Emission vector of wavelengths. \item ex Excitation vector of wavelengths.
#'   }
#'
#' @examples
#' library(eemR)
#' data("absorbance")
#'
#' folder <- system.file("extdata/cary/scans_day_1", package = "eemR")
#' eems <- eem_read(folder)
#' eems <- eem_extract(eems, "nano", remove = TRUE) # Remove the blank sample
#'
#' ## Remove scattering (1st order)
#' eems <- eem_remove_scattering(eems, "rayleigh")
#'
#' eems_corrected <- eem_inner_filter_effect(eems, absorbance = absorbance, pathlength = 1)
#'
#' op <- par(mfrow = c(2, 1))
#' plot(eems, which = 1)
#' plot(eems_corrected, which = 1)
#' par(op)
#'
#' @export
eem_inner_filter_effect <- function(eem, absorbance, pathlength = 1) {

  stopifnot(.is_eemlist(eem) | .is_eem(eem),
            is.data.frame(absorbance),
            is.numeric(pathlength))

  ## It is a list of eems, then call lapply
  if(.is_eemlist(eem)){

    res <- lapply(eem, eem_inner_filter_effect, absorbance = absorbance)

    class(res) <- class(eem)

    return(res)
  }

  #---------------------------------------------------------------------
  # Some checks
  #---------------------------------------------------------------------
  names(absorbance) <- tolower(names(absorbance))

  if(!any(names(absorbance) == "wavelength")){

    stop("'wavelength' variable was not found in the data frame.",
         call. = FALSE)
  }

  wl <- absorbance[, "wavelength"]

  if(!all(is_between(range(eem$em), min(wl), max(wl)))){

    stop("absorbance wavelenghts are not in the range of
         emission wavelengths", call. = FALSE)

  }

  if(!all(is_between(range(eem$ex), min(wl), max(wl)))){

    stop("absorbance wavelenghts are not in the range of
         excitation wavelengths", call. = FALSE)
  }

  spectra <- absorbance[, which(names(absorbance) == eem$sample)]

  ## absorbance spectra not found, we return the uncorected eem
  if(length(spectra) == 0){

    warning("Absorbance spectrum for ", eem$sample, " was not found. Returning uncorrected EEM.",
            call. = FALSE)

    return(eem)
  }

  #---------------------------------------------------------------------
  # Create the ife matrix
  #---------------------------------------------------------------------

  sf <- stats::splinefun(wl, spectra)

  ex <- sf(eem$ex)
  em <- sf(eem$em)

  total_absorbance <- sapply(ex, function(x){x + em})

  if(max(total_absorbance) > 1.5){
    cat("Total absorbance is > 1.5 (Atotal = ", max(total_absorbance), ")\n",
        "A 2-fold dilution is recommended. See ?eem_inner_filter_effect.",
        sep = "")
  }

  ife_correction_factor <- 10 ^ (-pathlength / 2 * (total_absorbance))

  cat("Range of IFE correction factors:",
      round(range(ife_correction_factor), digits = 4), "\n")

  cat("Range of total absorbance (Atotal) :",
      round(range(total_absorbance), digits = 4), "\n")

  x <- eem$x / ife_correction_factor

  ## Construct an eem object.
  res <- eem(file = eem$sample,
             x = x,
             ex = eem$ex,
             em = eem$em)

  attributes(res) <- attributes(eem)
  attr(res, "is_ife_corrected") <- TRUE

  return(res)

}

is_between <- function(x, a, b) {
  x >= a & x <= b
}

