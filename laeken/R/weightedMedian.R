# ------------------------------------------
# Authors: Andreas Alfons and Matthias Templ
#          Vienna University of Technology
# ------------------------------------------

#' Weighted median
#' 
#' Compute the weighted median (Eurostat definition).
#' 
#' The implementation strictly follows the Eurostat definition.
#' 
#' @param x a numeric vector.
#' @param weights an optional numeric vector giving the sample weights.
#' @param sorted a logical indicating whether the observations in \code{x} are
#' already sorted.
#' @param na.rm a logical indicating whether missing values in \code{x} should
#' be omitted.
#' @return The weighted median of values in \code{x} is returned.
#' 
#' @author Andreas Alfons and Matthias Templ
#' 
#' @seealso \code{\link{arpt}}, \code{\link{incMedian}},
#' \code{\link{weightedQuantile}}
#' 
#' @references Working group on Statistics on Income and Living Conditions
#' (2004) Common cross-sectional EU indicators based on EU-SILC; the gender pay
#' gap.  \emph{EU-SILC 131-rev/04}, Eurostat.
#' 
#' @keywords survey
#' 
#' @examples
#' data(eusilc)
#' weightedMedian(eusilc$eqIncome, eusilc$rb050)
#' 
#' @export

weightedMedian <- function(x, weights = NULL, sorted = FALSE, na.rm = FALSE) {
    weightedQuantile(x, weights, probs=0.5, sorted=sorted, na.rm=na.rm)
}
