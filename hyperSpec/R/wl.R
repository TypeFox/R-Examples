##' Getting and Setting the Wavelength Axis
##' \code{wl} returns the wavelength axis, \code{wl<-} sets it.
##' 
##' The wavelength axis of a \code{hyperSpec} object can be retrieved and
##' replaced with \code{wl} and \code{wl<-}, respectively.
##' 
##' When the wavelength axis is replaced, the colnames of \code{x@@data$spc} are
##' replaced by the rounded new wavelengths.  \code{digits} specifies the how
##' many significant digits should be used.
##' 
##' There are two ways to set the label of the new wavelength axis, see the
##' examples.  If no label is given, a warning will be issued.
##' 
##' @aliases wl
##' @param x a \code{hyperSpec} object
##' @return a numeric vector
##' @note \code{wl<-} always sets the complete wavelength axis, without
##'   changing the columns of the spectra matrix. If you rather want to cut the
##'   spectral range, use \code{\link[hyperSpec:extractreplace]{[}}, for
##'   interpolation along the spectral axis see
##'   \code{\link[hyperSpec]{spc.loess}} and for spectral binning
##'   \code{\link[hyperSpec]{spc.bin}}.
##' @author C. Beleites
##' @export
##' @seealso \code{\link[base]{signif}}
##' 
##' cutting the spectral range: \code{\link[hyperSpec:extractreplace]{[}}
##' 
##' interpolation along the spectral axis: \code{\link[hyperSpec]{spc.loess}}
##' 
##' spectral binning: \code{\link[hyperSpec]{spc.bin}}
##' @examples
##' 
##' 	wl (laser)
##' 
wl <- function (x){
  chk.hy (x)
  validObject (x)
  
  x@wavelength
}

###-----------------------------------------------------------------------------
###
###  .wl
###
###
".wl<-" <- function (x, digits = 6, value){
  x@wavelength <- value
  colnames (x@data$spc) <- signif (value, digits)

  x
}

##' @rdname wl
##' @export "wl<-"
##' @aliases wl<-
##' @usage
##' wl (x, label=NULL, digits=6) <- value
##' 
##' @param value either a numeric containing the new wavelength vector, or a
##'   list with \code{value$wl} containing the new wavelength vector and
##'   \code{value$label} holding the corresponding \code{label}.
##' @param label The label for the new wavelength axis. See \link{initialize}
##'   for details.
##' @param digits handed to \code{\link[base]{signif}}. See details.
##' @return \code{hyperSpec} object
##' @examples
##' 	# convert from wavelength to frequency 
##' 	plot (laser)
##' 	wl (laser, "f / Hz") <- 2.998e8 * wl (laser) * 1e9
##' 	plot (laser)
##' 	
##' 	# convert from Raman shift to wavelength
##' 	# excitation was at 785 nm
##' 	plot (chondro [1])
##' 	wl (chondro) <- list (wl = 1e7 / (1e7/785 - wl (chondro)), label = expression (lambda / nm))
##' 	plot (chondro [1])
##' 
"wl<-" <- function (x, label = NULL, digits = 6, value){

  chk.hy (x)
  validObject (x)

  if (is.list (value)){
    label <- value$label
    value <- value$wl
  }

  .wl (x) <- value

  x@label$.wavelength <- label

  validObject (x)

  x
}
