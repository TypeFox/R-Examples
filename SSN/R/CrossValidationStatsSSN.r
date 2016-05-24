CrossValidationStatsSSN <- function(object)
{
	cv.out <- CrossValidationSSN(object)
	zt <- object$sampinf$z
	n <-  object$sampinf$obs.sample.size
	dfobs <- n - object$sampinf$rankX
	cv.stats <- data.frame(bias = sum(cv.out[,"cv.pred"] - zt)/n,
		std.bias = sum((cv.out[,"cv.pred"] - zt)/sqrt(cv.out[,"cv.se"]))/n,
		RMSPE = sqrt(sum((zt - cv.out[,"cv.pred"])^2)/n),
		RAV = sqrt(sum(cv.out[,"cv.se"]^2)/n),
		std.MSPE = sum((cv.out[,"cv.pred"] - zt)^2/cv.out[,"cv.se"]^2)/n,
		cov.80 = sum(abs((zt - cv.out[,"cv.pred"])/cv.out[,"cv.se"]) < qt(.90, dfobs))/n,
		cov.90 = sum(abs((zt - cv.out[,"cv.pred"])/cv.out[,"cv.se"]) < qt(.95, dfobs))/n,
		cov.95 = sum(abs((zt - cv.out[,"cv.pred"])/cv.out[,"cv.se"]) < qt(.975, dfobs))/n)
  cv.stats
}

