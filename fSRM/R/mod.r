#' @title Get modification indices for a fSRM object
#' @aliases mod
#'
#' @description
#' Get modification indices for a fSRM object.
#'
#' @export
#' @param x A fSRM object.
#' @param minMI Minimum size of modification indices to be printed.

mod <- function(x, minMI = 10) {
	if ((x$means == TRUE | x$diff == TRUE) & packageVersion("lavaan") <= "0.5.15") {
		stop("Modification indices do not work when mean structure or delta method are used. If you upgrade to lavaan 0.5.16 (or later), you can get MI for these models as well.")
	}
	
	MI <- modindices(x$fit, standardized=TRUE)
	MI <- MI[order(MI$mi, decreasing=TRUE), ]
	
	# Joereskog: MI > 5 before consideration of respecification
	# Rosseel: MI > 10 before consideration of respecification
	if (max(MI$mi, na.rm=TRUE) < minMI) {
		print(paste0("No modification index is larger then ", minMI, "."))
		invisible(NULL)
	}
	MI2 <- MI[!is.na(MI$mi) & MI$mi>minMI, ]
	if (nrow(MI2)>0) {
		return(MI2)
	} else {
		invisible(NULL)
	}
}