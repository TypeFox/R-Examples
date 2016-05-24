##' Convert Hz to the mel scale
##' 
##' The calculation is done using the formulae mel = 1/log(2) * (log(1 +
##' (Hz/1000))) * 1000 where Hz is the frequency in Hz.
##' 
##' If 'data' is a spectral object, then the frequencies are changed so that
##' they are proportional to the mel scale and such that the mel intervals
##' between frequencies are constant between the lowest and highest
##' frequencies. More specifically, suppose that a spectral object has
##' frequencies at 0, 1000, 2000, 3000, 4000 Hz. Then the corresponding
##' frequencies extend in mel between 0 and 2321.928 mel (=4000 Hz in mels) in
##' four equal intervals, and linear interpolation is used with the 'approx'
##' function to obtain the dB values at those frequencies.
##' 
##' @aliases mel mel.trackdata mel.spectral
##' @param a A vector or matrix of data or a spectral object.
##' @return A vector or matrix or spectral object of the same length and
##' dimensions as data.
##' @author Jonathan Harrington
##' @seealso \code{\link{bark}}, \code{\link{plot.spectral}}
##' @references Traunmueller, H. (1990) \"Analytical expressions for the
##' tonotopic sensory scale\" J. Acoust. Soc. Am. 88: 97-100.
##' @keywords math
##' @examples
##' 
##' 
##'  #convert Hertz values to mel
##' 
##'  vec <- c(500, 1500, 2500)
##'  vec
##'  mel(vec)
##' 
##'         
##'         # convert Hertz values to mel 
##' 
##'         mel(vec)
##'         
##' 
##'         # convert the \$data values in a trackdata object to mel
##'         # create a new track data object
##' 
##'         t1 <- dip.fdat
##'         t1[1]
##'         
##'         # convert Hertz to mel
##' 
##'         t1$data <- mel(t1$data)
##'         t1[1]
##' 
##' # warp the frequency axis of a spectral object such
##' # that it is proportional to the mel scale.
##' 
##' w = mel(e.dft)
##' par(mfrow=c(1,2))
##' plot(w, type="l")
##' 
##' 
##' # The values of w are at equal mel intervals. Compare
##' # with
##' 
##' plot(e.dft, freq=mel(trackfreq(e.dft)))
##' 
##' # the latter has a greater concentration of values
##' # in a higher frequency range.
##' 
##' @export mel
"mel" <- function(a)
{
  UseMethod("mel")
}









##' mel default
##' 
##' see function
##' 
##' 
##' @keywords internal
##' @export
"mel.default" <- function (a) 
{
  1/log(2) * (log(1 + (a/1000))) * 1000
}
