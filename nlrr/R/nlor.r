#' Odds ratio plot for dose - response non-linear continuous exposure.
#'
#' Calculates non-linear odds ratio and plot OR vs. a continuous variable.
#'
#' @param outcome the outcome variable
#' @param exposure the exposure variable
#' @param covar a covariats list
#' @param ref reference value for the continuous variable
#' @param knum number of knots
#' @param data name of a dataset
#'
#' @export
#' @import rms
#' @import Hmisc
#' @importFrom stats binomial coef glm na.exclude vcov
#' @examples
#' sum1 <- nlor('dm', 'lipid', covar = c('age', 'gender'), 0.6, data = Lipid)
#' head(sum1)


nlor <- function(outcome, exposure, covar = NULL, ref = NULL, knum = 4, data){

	#data checking
	colname <- colnames(data)

	colnum <- which(colname %in% c(outcome, exposure, covar))

	data <- na.exclude(data[, colnum])

	if (length(covar) == 0) {

		covarM = matrix(0, nrow = dim(data)[1], ncol = 1)

	}
	else {

		covarM <- as.matrix(data[, covar])

	}

	#restricted cubic spline
	xs <- rcspline.eval(data[, exposure], nk = knum, inclx = TRUE)

	fit <- glm(data[, outcome] ~ xs + covarM, family = binomial)

	beta1 <- coef(summary(fit))[2:knum,1]

	xs <- xs[order(xs[,1]),]

	#reference value
	if (length(ref) == 0) {

		base0 <- xs[round(dim(xs)[1]/2),]

	}
	else {

		base0 <- xs[which.min(abs(xs[,1] - ref)),]

	}

	xref <- base0[1]

	#difference from reference value
	diff2ref <- t(t(xs) - base0)

	# log odds ratio
	logor <- diff2ref %*% beta1

	# se of logor
	varcov <- vcov(fit)[2:knum,2:knum]

	selogor <- sqrt(diag(diff2ref %*% varcov %*% t(diff2ref)))

	#Odds Ratio and 95%CI
	or <- exp(logor)

	orup <- exp(logor + 1.96 * selogor)

	orlow <- exp(logor - 1.96 * selogor)

	newdata <- cbind(exposure = xs[,1], or = or, orse = selogor, orup = orup, orlow = orlow, xref = xref)

	colnames(newdata) <- c(exposure, 'or', 'orse', 'orup', 'orlow', 'xref')

	return(newdata)

}




