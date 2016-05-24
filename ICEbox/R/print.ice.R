print.ice = function(x, ...){
	cat("ice object generated on data with n = ", nrow(x$ice_curves), " for predictor \"", x$predictor, "\"\n", sep = "")
	cat("predictor considered ", ifelse(x$nominal_axis, "discrete", "continuous"), ", logodds ", ifelse(x$logodds, "on", "off"), "\n", sep = "")
}
