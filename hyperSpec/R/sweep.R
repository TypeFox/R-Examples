.sweep <- function (x, MARGIN, STATS, FUN = "-",
										check.margin = TRUE, ...){
	validObject (x)
	
	if (is (STATS, "hyperSpec")){
		validObject (STATS)
		STATS <- STATS@data$spc
	} else if (is (STATS, "function")) {
		STATS <- apply (x, MARGIN, STATS)@data$spc
	}
	
	x@data$spc <- sweep (x = x@data$spc, MARGIN = MARGIN, STATS = STATS,
											 FUN = FUN, check.margin = check.margin, ...)
	
	x
}

##' Sweep Summary Statistic out of an hyperSpec Object
##' \code{\link[base]{sweep}} for \code{hyperSpec} objects.
##' 
##' Calls \code{\link[base]{sweep}} for the spectra matrix.
##' 
##' \code{sweep} is useful for some spectra preprocessing, like offset
##' correction, substraction of background spectra, and normalization of the
##' spectra.
##'
##' @name sweep
##' @rdname sweep
##' @aliases sweep-methods sweep,hyperSpec-method
##' @docType methods
##' @param x a \code{hyperSpec object.}
##' @param MARGIN direction of the spectra matrix that \code{STATS} goees
##'   along.
##' @param STATS the summary statistic to sweep out. Either a vector or a
##'   \code{hyperSpec} object.
##' 
##' hyperSpec offers a non-standard convenience function: if \code{STATS} is a
##'   function, this function is applied first (with the same \code{MARGIN}) to
##'   compute the statistic. However, no further arguments to the apply
##'   function can be given.  See the examples.
##' @param FUN the function to do the sweeping, e.g. `-` or `/`.
##' @param check.margin If \code{TRUE} (the default), warn if the length or
##'   dimensions of \code{STATS} do not match the specified dimensions of
##'   \code{x}.  Set to \code{FALSE} for a small speed gain when you
##'   \emph{know} that dimensions match.
##' @param ... further arguments for \code{FUN}
##' @return A \code{hyperSpec} object.
##' @author C. Beleites
##' @seealso \code{\link[base]{sweep}}
##' @keywords methods
##' @export 
##' @examples
##' 
##' ## Substract the background / slide / blank spectrum
##' # the example data does not have spectra of the empty slide,
##' # so instead the overall composition of the sample is substracted
##' background <- apply (chondro, 2, quantile, probs = 0.05)
##' corrected <- sweep (chondro, 2, background, "-")
##' plot (corrected, "spcprctl5")
##' 
##' ## Offset correction
##' offsets <- apply (chondro, 1, min)
##' corrected <- sweep (chondro, 1, offsets, "-")
##' plot (corrected, "spcprctl5")
##' 
##' ## Min-max normalization (on max amide I)
##' # the minimum is set to zero by the offset correction.
##' factor <- apply (corrected, 1, max)
##' mm.corrected <- sweep (corrected, 1, factor, "/")
##' plot (mm.corrected, "spcprctl5")
##' 
##' ## convenience: give function to compute STATS:
##' mm.corrected2 <- sweep (corrected, 1, max, "/")
##' plot (mm.corrected2)
##' 
##' ## checking
##' stopifnot (all (mm.corrected2 == mm.corrected))
##' 
setMethod ("sweep", signature = signature (x = "hyperSpec"), .sweep)
