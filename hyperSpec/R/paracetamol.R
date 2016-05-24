##' Paracetamol Spectrum
##' A Raman spectrum of a paracetamol tablet.
##' 
##' 
##' @name paracetamol
##' @docType data
##' @format The spectrum was acquired with a Renishaw InVia spectrometer from
##'   100 to 3200 cm^-1 in step scan mode. Thus the spectrum has several
##'   overlapping wavelength regions.
##' @author C. Beleites
##' @keywords datasets
##' @examples
##' 
##' paracetamol
##' 
##' plot (paracetamol)
##' plotspc (paracetamol, c (min ~ 1750, 2800 ~ max), xoffset = 800,
##' wl.reverse = TRUE)
##' 
NULL

