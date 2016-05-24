##' Function that applies a function to an EMU spectral object
##' 
##' Applies a function to an EMU spectral object.
##' 
##' fapply performs a similar operation to apply except that it is specifically
##' designed for handling EMU spectral objects.
##' 
##' @param specdata A matrix or trackdata object of class spectral
##' @param fun A function to be applied.
##' @param \dots Optional arguments to fun
##' @param power A single element logical vector. If T, convert specdata to
##' power values i.e. apply the function to a * specdata
##' \eqn{\mbox{\textasciicircum}}{^}b or a * specdata\$data
##' \eqn{\mbox{\textasciicircum}}{^}b where a and b powcoeffs defined below.
##' @param powcoeffs A 2 element numeric vector for converting dB values to
##' power values. Defaults to a = 10 and b = 10. See \code{power}.
##' @return If the output has the same dimensions has the input, then an object
##' of the same dimensionality and class is returned. Otherwise it may be a
##' vector or matrix depending on the function that is applied.  ...
##' @section Warning : The function can be very slow if applied to a large
##' trackdata object. In this case, it may be faster to use a for-loop with the
##' desired function around \$data
##' @author Jonathan Harrington
##' @seealso \code{\link{apply}} \code{\link{by.trackdata}}
##' @keywords utilities
##' @examples
##' 
##' # mean value per spectrum, input is a spectral matrix
##' m <- fapply(vowlax.dft.5, sapply, FUN=mean)
##' # as above but after converting dB to powers before
##' # applying the function
##' m <- fapply(vowlax.dft.5, sapply, FUN=mean, power=TRUE)
##' # spectral range
##' r <- fapply(vowlax.dft.5, range)
##' # spectral moments applied to a trackdata object
##' # m is a four-dimensional trackdata object
##' m <- fapply(fric.dft, moments)
##' # 1st 3 DCT coefficients calculated in a spectral matrix
##' # d is a 3-columned matrix
##' d <- fapply(vowlax.dft.5, dct, 3)
##' # dct-smooth with 10 coefficients. d2 is spectral matrix
##' d2 <- fapply(vowlax.dft.5, dct, 10, TRUE)
##' # dct-smooth a trackdata object with 10 coefficients
##' d3 <- fapply(fric.dft[1:4,], dct, 10, TRUE)
##' 
##' 
##' @export fapply
'fapply' <- function (specdata, fun, ..., power = FALSE, 
                      powcoeffs = c(10, 10)) 
{
  if (!is.spectral(specdata)) 
    stop("object must be of class spectral")
  if (power) 
    specdata <- dbtopower(specdata, powcoeffs[1], powcoeffs[2])
  dnames <- dimnames(specdata)
  omat <- NULL
  if (is.trackdata(specdata)) {
    class(specdata$data) <- "spectral"
    for (j in 1:nrow(specdata$data)) {
      vals <- fun(specdata$data[j, ], ...)
      omat <- rbind(omat, vals)
    }
    if (ncol(omat) == ncol(specdata)) {
      dimnames(omat) <- dnames
      if (!is.spectral(omat)) 
        omat <- as.spectral(omat, attr(specdata$data, 
                                       "fs"))
    }
    else dimnames(omat) <- list(dnames[[1]], NULL)
    if (power) 
      omat <- dbtopower(omat, powcoeffs[1], powcoeffs[2], 
                        inv = TRUE)
    return(as.trackdata(omat, specdata$index, specdata$ftime))
  }
  else {
    if(!is.matrix(specdata))
    {
      samfreq = max(trackfreq(specdata))*2
      specdata = as.spectral(rbind(specdata), samfreq)
    }
    
    for (j in 1:nrow(specdata)) {
      vals <- fun(specdata[j, ], ...)
      omat <- rbind(omat, vals)
    }
    if (ncol(omat) == ncol(specdata)) {
      dimnames(omat) <- dnames
      if (!is.spectral(omat)) {
        class(omat) <- c(class(omat), "spectral")
        attr(omat, "fs") <- attr(specdata, "fs")
      }
    }
    else dimnames(omat) <- list(dnames[[1]], NULL)
    if (power) 
      omat <- dbtopower(omat, powcoeffs[1], powcoeffs[2], 
                        inv = TRUE)
    return(omat)
  }
}
