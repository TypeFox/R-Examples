#' Odds ratio plot for dose - response non-linear continuous exposure.
#'
#' Calculates non-linear odds ratio and plot OR vs. a continuous variable.
#'
#' @param exposure the exposure variable
#' @param or odds ratio
#' @param xlab x-axis
#' @param data name of a dataset
#'
#' @export
#' @importFrom graphics abline axis lines plot polygon
#' @examples
#' sum1 <- nlor('dm', 'lipid', covar = c('age', 'gender'), 0.6, data = Lipid)
#' head(sum1)
#' nlorplot('lipid', 'or', data = sum1, xlab = 'Lipid')

nlorplot <- function(exposure, or, data, xlab = NULL){

	#y-axis range
	ymin <- min(data[, 'orlow'])
	ymax <- max(data[, 'orup'])

	#xlab
	if (length(xlab) == 0){

		xlab1 <- c(exposure)

	}
	else {

		xlab1 <- xlab

	}

	#plot
	plot(data[, exposure], data[, or], type = 'l',
			ylim = c(ymin, ymax), ylab = 'OR (95% CI)', xlab = xlab1)

	lines(data[, exposure], data[, 'orup'], col = 'grey')

	lines(data[, exposure], data[, 'orlow'], col = 'grey')

	polygon(c(data[, exposure], rev(data[, exposure])),
				c(data[, 'orup'], rev(data[, 'orlow'])), col = 'grey', border = 'grey')

	lines(data[, exposure], data[, or], col = 'black')

	abline(h = 1, v = data[1, 'xref'], lty = 2)

	axis(side = 1, at = data[1, 'xref'], labels = round(data[1, 'xref'], 2))

}
