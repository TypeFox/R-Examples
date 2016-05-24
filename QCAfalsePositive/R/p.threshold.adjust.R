
#' Calculate Vector of p-value Thresholds for Multiple Inference
#'
#' Internal function. Calculates the adjusted thresholds for multiple significance
#' tests, assuming that the original threshold for a single test is p<0.05.
#' Used by \code{\link{fsQCApermTest}} to calculate confidence intervals.
#' @param total.configs The total number of hypotheses tested, or the number of 
#' configurations utilized by the Quine-McCluskey algorithm in fsQCA (including 
#' logical remainders, if they are used in the analysis).
#' @param my.method The adjustment method used to calculate p-values (see 
#' \code{\link{p.adjust}} for details).
#' @return Numeric vector giving adjusted p-value thresholds, from smallest to
#' largest.
#' @keywords p-value threshold fsQCA
#' @export
#' @examples
#' p.threshold.adjust(10, "holm")


p.threshold.adjust <- function(total.configs, my.method){
	pvec <- NULL
	for(m in 1:total.configs){
		pvec <- c(pvec, p.adjust(c(rep(1, m-1), 0.05, rep(0, total.configs-m)), method=my.method)[m])
	}
	out.vec <- rep(0.05, total.configs) / (pvec/0.05)
	rev(out.vec)
}
