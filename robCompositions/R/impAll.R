#' Replacement of rounded zeros and missing values.
#' 
#' Parametric replacement of rounded zeros and missing values for compositional
#' data using classical and robust methods based on ilr-transformations with
#' special choice of balances. Values under detection limit should be saved
#' with the negative value of the detection limit (per variable). Missing
#' values should be coded as NA.
#' 
#' This is a wrapper function that calls \emph{impRZilr()} for the replacement
#' of zeros and \emph{impCoda} for the imputation of missing values
#' sequentially. The detection limit is automatically derived form negative
#' numbers in the data set.
#' 
#' @param x data frame
#' @return The imputed data set.
#' @note This function is mainly used by the compositionsGUI.
#' @author Jiri Eichler
#' @seealso \code{\link{impCoda}}, \code{\link{impRZilr}}
#' @export
#' @references Hron, K. and Templ, M. and Filzmoser, P. (2010) Imputation of
#' missing values for compositional data using classical and robust methods,
#' \emph{Computational Statistics and Data Analysis}, vol 54 (12), pages
#' 3095-3107.
#' 
#' Martin-Fernandez, J.A.  and Hron, K. and Templ, M. and Filzmoser, P. and
#' Palarea-Albaladejo, J. (2012) Model-based replacement of rounded zeros in
#' compositional data: Classical and robust approaches, \emph{Computational
#' Statistics}, 56 (2012), S. 2688 - 2704.
#' @keywords manip
#' @examples
#' 
#' ## see the compositionsGUI
#' 
impAll <-
function(x) {
	maxLimits = apply(x, 2, min, na.rm = TRUE)
	maxLimits = -maxLimits
	maxLimits[maxLimits < 0] = 0

	if(any(is.na(x))) {
		temp <- x
		for(i in 1:length(maxLimits)) {
			temp[na.omit(temp[i]) < 0, i] = maxLimits[i] * 2 / 3
		}

		res <- adjust(impCoda(temp))
		isna = is.na(x)
		x[isna] = res$xImp[isna]
	}
	if(any(x < 0)) {
		temp <- x
		temp[temp < 0] = 0

		res <- impRZilr(temp, dl=maxLimits)
		x = res$xImp
	}

	return(x)
}

