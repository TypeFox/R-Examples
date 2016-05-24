print.dice = function(x, ...){
	cat("dice object generated on data with n = ", nrow(x$d_ice_curves), " for predictor \"", x$predictor, "\"\n", sep = "")
	cat("predictor considered ", ifelse(x$nominal_axis, "discrete", "continuous"), ", logodds ", ifelse(x$logodds, "on", "off"), "\n", sep = "")
}
