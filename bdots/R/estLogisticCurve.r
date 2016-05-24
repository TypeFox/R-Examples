est.logistic.curve <- function(time, fixations, rho, params = NULL, cor = TRUE) {
	if(is.null(params)) {
		mini  <- find.mini(time,  fixations)
		peak  <- find.peak(time,  fixations)
		slope <- find.slope(time, fixations)
		cross <- find.cross(time, fixations)
	} else {
		mini  <- params[1]
		peak  <- params[2]
		slope <- params[3]
		cross <- params[4]
	}
	if(cor) {
		fit.curve <- tryCatch(gnls(fixations ~ mini + (peak - mini) / (1 + exp(4 * slope * (cross - (time)) / (peak - mini))), 
									start=c(mini = mini, peak = peak, slope = slope, cross = cross),
									correlation = corAR1(rho)), error = function(e) NULL)
		if(is.null(fit.curve)) cor <- FALSE
	}
	if(!cor) {
			fit.curve <- tryCatch(gnls(fixations ~ mini + (peak - mini) / (1 + exp(4 * slope * (cross - (time)) / (peak - mini))), 
									start=c(mini = mini, peak = peak, slope = slope, cross = cross)), error = function(e) NULL)
	}
	list(fit = fit.curve, cor = cor)
}