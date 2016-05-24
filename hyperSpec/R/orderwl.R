
##' Sorting the Wavelengths of a hyperSpec Object
##' Rearranges the \code{hyperSpec} object so that the wavelength vector is in increasing (or
##' decreasing) order.
##' 
##' The wavelength vector is sorted and the columns of the spectra matrix are rearranged accordingly.
##' 
##' @param x The \code{hyperSpec} object.
##' @param na.last,decreasing Handed to \code{\link[base]{order}}.
##' @return A \code{hyperSpec} object.
##' @author C. Beleites
##' @export
##' @seealso \code{\link[base]{order}}
##' @examples
##' 
##' ## Example 1: different drawing order in plotspc
##' spc <- new ("hyperSpec", spc = matrix (rnorm (5) + 1:5, ncol = 5))
##' spc <- cbind (spc, spc+.5)
##' 
##' plot (spc, "spc")
##' text (wl (spc), spc [[]], as.character (1:10))
##' spc <- orderwl (spc)
##' plot (spc, "spc")
##' text (wl (spc), spc [[]], as.character (1:10))
##' 
##' ## Example 2
##' spc <- new ("hyperSpec", spc = matrix (rnorm (5)*2 + 1:5, ncol = 5))
##' spc <- cbind (spc, spc)
##' 
##' plot (seq_len(nwl(spc)), spc[[]], type = "b")
##' spc[[]]
##' 
##' spc <- orderwl (spc)
##' lines (seq_len(nwl(spc)), spc[[]], type = "l", col = "red")
##' spc[[]]
##' 
orderwl <- function (x, na.last = TRUE, decreasing = FALSE){
  chk.hy (x)
  validObject (x)

  .orderwl (x)
}

.orderwl <- function (x, na.last = TRUE, decreasing = FALSE){
  ord <- order (x@wavelength, na.last = na.last, decreasing = decreasing)

  if (any (ord != seq_along (x@wavelength))){
    x@data$spc <-  x@data$spc [, ord, drop = FALSE]
    .wl(x) <- x@wavelength [ord]
  }

  x
}
