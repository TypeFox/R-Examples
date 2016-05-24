#---------------------------------------------------------------------
# Warning message when there is a mismatch between metrics wavelengths
# and data wavelengths.
#---------------------------------------------------------------------
msg_warning_wavelength <- function(){
  msg <- "This metric uses either excitation or emission wavelenghts that were not present in the data. Data has been interpolated to fit the requested wavelengths."
  return(msg)
}

#' Calculate the fluorescence index (FI)
#'
#' @template template_eem
#'
#' @template template_section_interp2
#'
#' @references \url{http://doi.wiley.com/10.4319/lo.2001.46.1.0038}
#'
#' @return A data frame containing fluorescence index (FI) for each eem.
#' @export
#' @examples
#' file <- system.file("extdata/cary/scans_day_1/", "sample1.csv", package = "eemR")
#' eem <- eem_read(file)
#'
#' eem_fluorescence_index(eem)

eem_fluorescence_index <- function(eem, verbose = TRUE){

  stopifnot(.is_eemlist(eem) | .is_eem(eem))

  ## It is a list of eems, then call lapply
  if(.is_eemlist(eem)){

    res <- lapply(eem, eem_fluorescence_index, verbose = verbose)
    res <- dplyr::bind_rows(res)

    return(res)
  }

  if(!all(370 %in% eem$ex & c(450, 500) %in% eem$em) & verbose){

    warning(msg_warning_wavelength(), call. = FALSE)

  }

  fluo_450 <- pracma::interp2(eem$ex, eem$em, eem$x, 370, 450)
  fluo_500 <- pracma::interp2(eem$ex, eem$em, eem$x, 370, 500)

  fi <- fluo_450 / fluo_500

  return(data.frame(sample = eem$sample, fi = fi, stringsAsFactors = FALSE))

}

#' Extract fluorescence peaks
#'
#' @template template_eem
#'
#' @template template_section_interp2
#'
#' @return A data frame containing peaks B, T, A, M and C for each eem. See
#'   details for more information.
#'
#' @details According to Coble (1996), peaks are defined as follow:
#'
#'   Peak B: ex = 275 nm, em = 310 nm
#'
#'   Peak T: ex = 275 nm, em = 340 nm
#'
#'   Peak A: ex = 260 nm, em = 380:460 nm
#'
#'   Peak M: ex = 312 nm, em = 380:420 nm
#'
#'   peak C: ex = 350 nm, em = 420:480 nm
#'
#'   Given that peaks A, M and C are not defined at fix emission wavelength,
#'   the maximum fluorescence value in the region is extracted.
#'
#' @references Coble, P. G. (1996). Characterization of marine and terrestrial
#'   DOM in seawater using excitation-emission matrix spectroscopy. Marine
#'   Chemistry, 51(4), 325-346.
#'
#'   \url{http://doi.org/10.1016/0304-4203(95)00062-3}
#'
#' @examples
#' file <- system.file("extdata/cary/scans_day_1/", "sample1.csv", package = "eemR")
#' eem <- eem_read(file)
#'
#' eem_coble_peaks(eem)
#'
#' @export
eem_coble_peaks <- function(eem, verbose = TRUE){

  stopifnot(.is_eemlist(eem) | .is_eem(eem))

  ## It is a list of eems, then call lapply
  if(.is_eemlist(eem)){

    res <- lapply(eem, eem_coble_peaks, verbose = verbose)
    res <- dplyr::bind_rows(res)

    return(res)
  }

  coble_ex_peak <- list(b = 275, t = 275, a = 260, m = 312, c = 350)

  if(!all(coble_ex_peak %in% eem$ex) & verbose){

    warning(msg_warning_wavelength(), call. = FALSE)

  }

  ## Get the peaks
  b <- pracma::interp2(eem$ex, eem$em, eem$x, 275, 310)

  t <- pracma::interp2(eem$ex, eem$em, eem$x, 275, 340)

  a <- max(pracma::interp2(eem$ex, eem$em, eem$x,
                           rep(260, length(380:460)), 380:460))

  m <- max(pracma::interp2(eem$ex, eem$em, eem$x,
                           rep(312, length(380:420)), 380:420))

  c <- max(pracma::interp2(eem$ex, eem$em, eem$x,
                           rep(350, length(420:480)), 420:480))

  #--------------------------------------------
  # Return the data
  #--------------------------------------------
  return(data.frame(sample = eem$sample,
                    b = b,
                    t = t,
                    a = a,
                    m = m,
                    c = c,
                    stringsAsFactors = FALSE))
}

#' Calculate the fluorescence humification index (HIX)
#'
#' @template template_eem
#' @param scale Logical indicating if HIX should be scaled, default is FALSE.
#'   See details for more information.
#'
#' @template template_section_interp2
#'
#' @description The fluorescence humification index (HIX), which compares two
#'   broad aromatic dominated fluorescence maxima, is calculated at 254 nm
#'   excitation by dividing the sum of fluorescence intensities between emission
#'   435 to 480 nm by the the sum of fluorescence intensities between 300 to 345
#'   nm.
#'
#' @references Ohno, T. (2002). Fluorescence Inner-Filtering Correction for
#'   Determining the Humification Index of Dissolved Organic Matter.
#'   Environmental Science & Technology, 36(4), 742-746.
#'
#'   \url{http://doi.org/10.1021/es0155276}
#'
#' @return A data frame containing the humification index (HIX) for each eem.
#' @export
#' @examples
#' file <- system.file("extdata/cary/scans_day_1/", package = "eemR")
#' eem <- eem_read(file)
#'
#' eem_humification_index(eem)
#'
eem_humification_index <- function(eem, scale = FALSE, verbose = TRUE) {

  stopifnot(.is_eemlist(eem) | .is_eem(eem),
            is.logical(scale))

  ## It is a list of eems, then call lapply
  if(.is_eemlist(eem)){

    res <- lapply(eem, eem_humification_index, verbose = verbose, scale = scale)
    res <- dplyr::bind_rows(res)

    return(res)
  }

  #---------------------------------------------------------------------
  # Get the data and calculate the humification index (HIX)
  #---------------------------------------------------------------------

  if(!254 %in% eem$ex & verbose){

    warning(msg_warning_wavelength(), call. = FALSE)

  }

  em_435_480 <- seq(from = 435, to = 480, by = 1)
  em_300_345 <- seq(from = 300, to = 345, by = 1)
  ex_254 <- rep(254, length(em_300_345))

  sum_em_435_480 <- sum(pracma::interp2(eem$ex, eem$em, eem$x,
                                        ex_254, em_435_480))

  sum_em_300_345 <- sum(pracma::interp2(eem$ex, eem$em, eem$x,
                                        ex_254, em_300_345))

  if(scale){

    hix <- sum_em_435_480 / (sum_em_300_345 + sum_em_435_480)

  }else{

    hix <- sum_em_435_480 / sum_em_300_345

  }

  return(data.frame(sample = eem$sample, hix = hix, stringsAsFactors = FALSE))
}

#' Calculate the biological fluorescence index (BIX)
#'
#' @template template_eem
#'
#' @template template_section_interp2
#'
#' @description The biological fluorescence index (BIX) is calculated by
#'   dividing the fluorescence at excitation 310 nm and emission at 380 nm (ex =
#'   310, em = 380) by that at excitation 310 nm and emission at 430 nm (ex =
#'   310, em = 430).
#'
#' @references Huguet, A., Vacher, L., Relexans, S., Saubusse, S., Froidefond,
#'   J. M., & Parlanti, E. (2009). Properties of fluorescent dissolved organic
#'   matter in the Gironde Estuary. Organic Geochemistry, 40(6), 706-719.
#'
#'   \url{http://doi.org/10.1016/j.orggeochem.2009.03.002}
#'
#' @return A data frame containing the biological index (BIX) for each eem.
#' @export
#' @examples
#' file <- system.file("extdata/cary/scans_day_1/", package = "eemR")
#' eem <- eem_read(file)
#'
#' eem_biological_index(eem)
#'
eem_biological_index <- function(eem, verbose = TRUE) {

  stopifnot(.is_eemlist(eem) | .is_eem(eem))

  ## It is a list of eems, then call lapply
  if(.is_eemlist(eem)){

    res <- lapply(eem, eem_biological_index, verbose = verbose)
    res <- dplyr::bind_rows(res)

    return(res)
  }

  #---------------------------------------------------------------------
  # Get the data and calculate the biological index (BIX)
  #---------------------------------------------------------------------

  if(!all(310 %in% eem$ex & c(380, 430) %in% eem$em) & verbose){
    warning(msg_warning_wavelength(), call. = FALSE)
  }

  fluo_380 <- pracma::interp2(eem$ex, eem$em, eem$x, 310, 380)
  fluo_430 <- pracma::interp2(eem$ex, eem$em, eem$x, 310, 430)

  bix <- fluo_380 / fluo_430

  return(data.frame(sample = eem$sample, bix = bix, stringsAsFactors = FALSE))
}
