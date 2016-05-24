##' Convert Hertz to Bark and Bark to Hertz
##' 
##' The calculation is done using the formulae Traunmueller (1990)
##' 
##' If 'data' is a spectral object, then
##' 
##' the frequencies are changed so that they are proportional
##' 
##' to the Bark scale and such that the Bark intervals
##' 
##' between frequencies are con stant between the lowest
##' 
##' and highest frequencies. More specifically,
##' 
##' suppose that a spectral object has frequencies
##' 
##' at 0, 1000, 2000, 3000, 4000 Hz. Then the corresponding
##' 
##' frequencies extend in Bark between 0 and 17.46329 Bark
##' 
##' in four equal intervals, and linear interpolation
##' 
##' is used with the 'approx' function to obtain
##' 
##' the dB values at those frequencies. Negative frequencies
##' 
##' which are obtained for values of about less than 40 Hz
##' 
##' are removed in the case of spectral objects.
##' 
##' @aliases bark bark.trackdata bark.spectral bark.default
##' @param f A vector or matrix of data or a spectral object.
##' @param inv A single element logical vector. If F, data are converted from
##' Hertz to Bark, if T, data are converted from Bark to Hertz. (Does not apply
##' if 'data' is an oject of class 'spectral'.
##' @param \dots for generic only
##' @return
##' 
##' A vector or matrix or spectral object of the same length and dimensions as
##' data.
##' @author Jonathan Harrington
##' @seealso
##' 
##' \code{\link{mel}},
##' 
##' \code{\link{plot.spectral}}
##' @references Traunmueller, H. (1990) "Analytical expressions for the
##' tonotopic sensory scale" J. Acoust. Soc. Am. 88: 97-100.
##' @keywords math
##' @examples
##' 
##' 
##' 
##' 
##'    # convert Hertz values to Bark
##' 
##'         vec <- c(500, 1500, 2500)
##' 
##'         vec
##' 
##'         bark(vec)
##' 
##' 
##' 
##'         
##' 
##'         # convert Hertz values to Bark and back to Hertz
##' 
##'         bark(bark(vec, inv=TRUE))
##' 
##'         
##' 
##' 
##' 
##'         # convert the \$data values in a trackdata object to Bark
##' 
##'         # create a new track data object
##' 
##'         t1 <- dip.fdat
##' 
##'         t1[1]
##' 
##'         
##' 
##' 
##' 
##'         # convert Hertz to Bark
##' 
##'         t1$data <- bark(t1$data)
##' 
##'         t1[1]
##' 
##' 
##' 
##' # warp the frequency axis of a spectral object such
##' 
##' # that it is proportional to the Bark scale.
##' 
##' w = bark(e.dft)
##' 
##' par(mfrow=c(1,2))
##' 
##' plot(w, type="l")
##' 
##' 
##' 
##' # The values of w are at equal Bark intervals. Compare
##' 
##' # with
##' 
##' plot(e.dft, freq=bark(trackfreq(e.dft)))
##' 
##' # the latter has a greater concentration of values
##' 
##' # in a higher frequency range.
##' 
##' 
##' 
##' 
##' @export bark
"bark" <- function (f, inv = FALSE, ...) {
  UseMethod("bark") 
}


##' @export
"bark.default" <- function (f, inv = FALSE, ...) 
{
  if (!inv) {
    result = ((26.81 * f)/(1960 + f)) - 0.53
  } else {
    result = (1960 * (f + 0.53))/(26.28 - f)
  }
  return(result)
}
