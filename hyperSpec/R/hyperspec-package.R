##' Interface for hyperspectral data sets
##' This package gives an interface to handle hyperspectral data sets in R.
##' Hyperspectral data are spatially or time-resolved spectra, or spectra with
##' any other kind of information associated with the spectra. E.g. spectral
##' maps or images, time series, calibration series, etc.
##' 
##' The spectra can be data as obtained in XRF, UV/VIS, Fluorescence, AES, NIR,
##' IR, Raman, NMR, MS, etc.
##' 
##' More generally, any data that is recorded over a discretized variable, e.g.
##' absorbance = f (wavelength), stored as a vector of absorbance values for
##' discrete wavelengths is suitable.
##' 
##' @name hyperSpec-package
##' @title Package hyperSpec
##' @docType package
##' @author C. Beleites
##' 
##' Maintainer: Claudia Beleites <chemometrie@@beleites.de>
##' @seealso \code{citation ("hyperSpec")} produces the correct citation.
##' 
##' \code{package?hyperSpec} for information about the package
##' 
##' \code{class?hyperSpec} for details on the S4 class provided by this
##'   package.
##' @rdname hyperSpec-package
##' @include flu.R
##' @include chondro.R
##' @include laser.R
##' @include paracetamol.R
##' @include barbiturates.R
##' @keywords package
if (!requireNamespace ("svUnit", quietly = TRUE)){
  `.test<-` <- function (f, value) {
      class (value) <-  c ("svTest", "function")
    attr (f, "test") <- value
    f
  }
} else {
 `.test<-` <- `test<-`
}

