#' Waveband constructor method
#'
#' Constructor for "waveband" objects that can be used as imput when calculating
#' irradiances.
#'
#' @param x any R object on which applying the function range yields an array of
#'   two numeric values, describimg a range of wavelengths (nm)
#' @param weight a character string "SWF" or "BSWF", use NULL (the defalt) to
#'   indicate no weighting used when calculating irradiance
#' @param SWF.e.fun a function giving multipliers for a spectral weighting
#'   function (energy) as a function of wavelength (nm)
#' @param SWF.q.fun a function giving multipliers for a spectral weighting
#'   function (quantum) as a function of wavelength (nm)
#' @param SWF.norm a numeric value giving the native normalization wavelength
#'   (nm) used by SWF.e.fun and SWF.q.fun
#' @param norm a single numeric value indicating the wavelength at which the SWF
#'   should be normalized to 1.0, in nm. "NULL" means no normalization.
#' @param hinges a numeric array giving the wavelengths at which the s.irrad
#'   should be inserted by interpolation, no interpolation is indicated by an
#'   empty array (numeric(0)), if NULL then interpolation will take place at
#'   both ends of the band.
#' @param wb.name character string giving the name for the waveband defined,
#'   default is NULL
#' @param wb.label character string giving the label of the waveband to be used
#'   for ploting, default is wb.name
#'
#' @return a \code{waveband} object
#'
#' @export
#' @examples
#' waveband(c(400,700))
#'
#' @family waveband constructors
#'
waveband <- function(x,
                     weight=NULL, SWF.e.fun=NULL, SWF.q.fun=NULL, norm=NULL,
                     SWF.norm=NULL, hinges=NULL, wb.name=NULL, wb.label=wb.name) {
  if (is.generic_spct(x) && is.null(wb.name)) {
    wb.name = "Total"
  }
  x.range <- range(x)
  new_waveband(x.range[1], x.range[2], weight=weight, SWF.e.fun=SWF.e.fun, SWF.q.fun=SWF.q.fun,
               norm=norm, SWF.norm=SWF.norm, hinges=hinges, wb.name=wb.name, wb.label=wb.label)
}

#' @describeIn waveband A less flexible variant
#' @param w.low numeric value, wavelength at the short end of the band (nm)
#' @param w.high numeric value, wavelength at the long end of the band (nm)
#'
#' @export
#' @exportClass waveband
#' @examples
#' new_waveband(400,700)
#'
new_waveband <- function(w.low, w.high,
                         weight=NULL, SWF.e.fun=NULL, SWF.q.fun=NULL, norm=NULL,
                         SWF.norm=NULL, hinges=NULL, wb.name=NULL, wb.label=wb.name){
  # we make sure that hinges is not NULL, as this would cause problems elsewhere
  # if we are not using a SWF then we do not need to add hinges as we will be anyway interpolating
  # raw irradiances rather than weighted irradiances
  if (is.null(hinges)) {
    hinges <- c(w.low - 1e-12, w.low, w.high - 1e-12, w.high)
  }
  if (!is.null(weight) && weight!="none") {
    #
    if (!is.null(SWF.e.fun) && is.null(SWF.q.fun)){
      if (!is.null(SWF.norm)){
        SWF.q.fun <- function(w.length){SWF.e.fun(w.length) *  SWF.norm / w.length}
      } else {
        warning("Warning: either both photon and energy SWFs should be supplied, or a value for the
              wavelength at which the function supplied is normalized should be supplied through SWF.norm")
      }
    } else if (!is.null(SWF.q.fun) && is.null(SWF.e.fun)  && !is.null(SWF.norm)){
      if (!is.null(SWF.norm)){
        SWF.e.fun <- function(w.length){SWF.q.fun(w.length) * w.length / SWF.norm}
      } else {
        warning("Warning: either both photon and energy SWFs should be supplied, or a value for the
              wavelength at which the function supplied is normalized should be supplied through SWF.norm")
      }
    } else if (is.null(SWF.e.fun) && is.null(SWF.q.fun)){
      warning("weight != NULL, but no SWFs supplied")
      return(NA)
    }
    if (is.null(wb.name)) {
      wb.name <- paste("range", as.character(signif(w.low, 4)), as.character(signif(w.high, 4)), "wtd", sep=".")
      wb.label <- wb.name
    }
  } else {
    weight <- "none"
    if (is.null(wb.name)) {
      wb.name <- paste("range", as.character(signif(w.low, 4)), as.character(signif(w.high, 4)), sep=".")
      wb.label <- wb.name
    }
  }
  w_band <- list(low = w.low, high = w.high,
                 weight = weight, SWF.e.fun = SWF.e.fun, SWF.q.fun = SWF.q.fun, SWF.norm = SWF.norm,
                 norm = norm, hinges = hinges, name = wb.name, label = wb.label)
  class(w_band) <- c("waveband", class(w_band))
  return(w_band)
}

#' List-of-wavebands constructor
#'
#' Build a list of unweighted "waveband" objects that can be used as imput when
#' calculating irradiances.
#'
#' @param x a numeric array of wavelengths to split at (nm), or a range of
#'   wavelengths or a generic_spct or a waveband.
#' @param list.names character vector with names for the component wavebands in
#'   the returned list (in order of increasing wavelength)
#' @param short.names logical indicating whether to use short or long names for
#'   wavebands
#' @param length.out numeric giving the number of regions to split the range
#'   into (ignored if w.length is not numeric).
#'
#' @return an un-named list of wabeband objects
#'
#' @export
#' @examples
#' split_bands(c(400,500,600))
#' split_bands(list(c(400,500),c(550,650)))
#' split_bands(list(A=c(400,500),B=c(550,650)))
#' split_bands(c(400,500,600), short.names=FALSE)
#' split_bands(c(400,500,600), list.names=c("a","b"))
#' split_bands(c(400,700), length.out=6)
#' split_bands(400:700, length.out=3)
#' split_bands(sun.spct, length.out=10)
#' split_bands(waveband(c(400,700)), length.out=5)
#'
#' @note \code{list.names} is used to assign names to the elements of the list,
#'   while the waveband objects themselves always retain their \code{wb.label}
#'   and \code{wb.name} as generated during their creation.
#'
#' @family waveband constructors
#'

split_bands <- function(x, list.names=NULL, short.names=is.null(list.names), length.out=NULL) {
  if (!is.any_spct(x) && !is.waveband(x) && is.list(x)) {
    x.len <- length(x)
    names.len <- length(list.names)
    if (names.len < x.len) {
      list.names <- names(x)
      names.len <- length(list.names)
    }
    if (names.len < x.len) {
      list.names <- paste("wb", letters[1:x.len], sep=".")
    }
    bands.out <- list()
    for (i in 1:x.len) {
      wb.temp <- split_bands(x[[i]],
                             list.names[i],
                             length.out = length.out)
      bands.out <- c(bands.out, wb.temp)
    }
    return(bands.out)
  }
  if (is.generic_spct(x) || is.waveband(x)) {
    w.length <- range(x)
  } else if (is.numeric(x)) {
    x <- unique(sort(x))
    if (length(x) > 1L) {
      w.length <- x
    } else {
      warning("At least two wavelength values are needed")
      return(list())
    }
  }
  wl.len <- length(w.length)
  if (!is.null(length.out) && is.numeric(length.out)) {
    if (length.out < 1L) {
      warning("'length.out' if non null must be >= 1")
      return(list())
    } else {
      w.length <- range(w.length)
      wl.len <- length.out + 1
      w.length <- seq(w.length[1], w.length[2], length.out=wl.len)
    }
  } else {
    length.out <- wl.len - 1
  }
  use.wb.names <- is.null(list.names)
  bands.out <- list()
  names <- character()
  for (i in 1:(wl.len - 1)) {
    wb.temp <- new_waveband(w.length[i], w.length[i+1],
                            hinges=NULL, wb.name=NULL)
    names[i] <- ifelse(use.wb.names, labels(wb.temp)[[1]], list.names[i])
    bands.out <- c(bands.out, list(wb.temp))
  }
  num.bands <- wl.len - 1

  if (length(bands.out) == 1 && is.null(list.names)) {
    return(bands.out[[1]])
  } else {
    if (short.names) {
      names(bands.out) <- paste("wb", 1:length(bands.out), sep="")
    } else {
      if (!is.null(list.names) && length(list.names) >= num.bands)  {
        names(bands.out) <- list.names[1:num.bands]
      } else {
        names(bands.out) <- names
      }
    }
    return(bands.out)
  }
}

#' Query if it is a waveband
#'
#' Functions to check if an object is waveband.
#'
#' @param x any R object
#'
#' @return is.waveband returns TRUE if its argument is a waveband and FALSE otherwise.
#'
#' @export
#'
is.waveband <- function(x) {
  inherits(x, "waveband")
}

### I need to add a check.waveband() method and use it in the constructor and maybe also
### add non-functional replacement operators.
###

check.waveband <- function(x, byref = FALSE, strict.range = FALSE, ...) {
  stopifnot(x[["low"]] < x[["high"]])
  stopifnot(x[["weight"]] == "none" && !(is.null(x[["SWF.e.fun"]] && is.null(x[["SWF.q.fun"]]))))
  stopifnot(x[["weight"]] != "none" && (is.null(x[["SWF.e.fun"]] || is.null(x[["SWF.q.fun"]]))))
  x
}
