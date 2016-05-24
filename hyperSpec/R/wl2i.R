###-----------------------------------------------------------------------------
###
### .getindex
###
###
## does the acual work of looking up the index for wl2i, .extract and .replace
## extrapolate = TRUE returns first resp. last index for wavelength outside hyperSpec@wavelength.
## extrapolate = FALSE returns NA in this case

.getindex <- function (x, wavelength, extrapolate = TRUE){
    if (! extrapolate) {
        wavelength [wavelength < min (x@wavelength)] <- NA
        wavelength [wavelength > max (x@wavelength)] <- NA
    }
    tmp <- wavelength [! is.na (wavelength)]
    if (length (tmp) > 0) {
        tmp <- sapply (tmp,
                         function (x, y) which.min (abs (x  - y)),
                         x@wavelength)
        wavelength [! is.na (wavelength)] <- tmp
    }
    wavelength
}

##' Conversion between Wavelength and Spectra Matrix Column
##' Index \code{wl2i} returns the column indices for the spectra matrix for the given wavelengths.
##' \code{i2wl} converts column indices into wavelengths.
##'
##' If \code{wavelength} is numeric, each of its elements is converted to the respective index.
##' Values outside the range of \code{x@@wavelength} become \code{NA}.
##' 
##' If the range is given as a formula (i.e. \code{start ~ end}, a sequence
##' 
##' index corresponding to start : index corresponding to end
##' 
##' is returned. If the wavelengths are not ordered, that may lead to chaos. In this case, call
##' \code{\link[hyperSpec]{orderwl}} first.
##' 
##' Two special variables can be used: \code{min} and \code{max}, corresponding to the lowest and
##' highest wavelength of \code{x}, respectively.
##' 
##' start and end may be complex numbers. The resulting index for a complex x is then
##' 
##' index (Re (x)) + Im (x)
##' 
##' @aliases wl2i
##' @param x a \code{hyperSpec} object
##' @param wavelength the wavelengths to be converted into column indices,
##'   either numeric or a formula, see details.
##' @param i the column indices into the spectra matrix for which the
##'   wavelength is to be computed
##' @return A numeric containing the resulting indices for \code{wl2i}
##' @author C. Beleites
##' @export
##' @examples
##' 
##' flu
##' wl2i (flu, 405 : 407)
##' wl2i (flu, 405 ~ 407)
##'   
##' ## beginning of the spectrum to 407 nm 
##' wl2i (flu, min ~ 407)
##' 
##' ## 2 data points from the beginning of the spectrum to 407 nm 
##' wl2i (flu, min + 2i ~ 407)
##' 
##' ## the first 3 data points   
##' wl2i (flu, min ~ min + 2i)
##' 
##' ## from 490 nm to end of the spectrum 
##' wl2i (flu, 490 ~ max)
##'   
##' ## the last 8 data points 
##' wl2i (flu, max - 7i ~ max)
##'   
##' ## get 450 nm +- 3 data points
##' wl2i (flu, 450 - 3i ~ 450 + 3i) 
##'   
##' wl2i (flu, 300 : 400) ## all NA: 
##' wl2i (flu, 600 ~ 700) ## NULL: completely outside flu's wavelength range
##' 
wl2i <- function (x, wavelength = stop ("wavelengths are required.")){
  chk.hy (x)
  validObject (x)

  ## special in formula
  max <- max (x@wavelength)
  min <- min (x@wavelength)

  envir <- attr (wavelength, ".Environment")

  `~` <- function (e1, e2){
    if (missing (e2))              # happens with formula ( ~ end)
      stop ("wavelength must be a both-sided formula")

    if (    (Re (e1) < min (x@wavelength) && Re (e2) < min (x@wavelength)) ||
        (Re (e1) > max (x@wavelength) && Re (e2) > max (x@wavelength))){
      NULL                       # wavelengths completely outside the wl. range of x
    } else {
      e1 <- .getindex (x, Re (e1)) + Im (e1)
      e2 <- .getindex (x, Re (e2)) + Im (e2)

      if (e1 <= 0 || e2 <= 0|| e1 > length (x@wavelength) || e2 > length (x@wavelength))
        warning ("wl2i: formula yields indices outside the object.")

      seq (e1, e2)
    }
  }

  .conv.range <- function (range){
    if (is.numeric (range)){
      .getindex (x, range, extrapolate = FALSE)
    } else
    eval (range)
  }

  if (is.list (wavelength)) {
    unlist (lapply (wavelength, .conv.range))
  } else {
    .conv.range (wavelength)
  }
}


##' @rdname wl2i
##' @aliases i2wl
##' @return \code{i2wl} returns a numeric with the wavelengths
##' @export
##' @examples
##'   
##' i2wl (chondro, 17:20)
##' 
i2wl <- function (x, i){
  chk.hy (x)
  validObject (x)

  x@wavelength [i]
}

