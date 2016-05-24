fitBspline <- function(dataValues, continuousCovariates = NULL, indicators = NULL, group = NULL, covariate = NULL) {
	parametersEstimated = sum(length(continuousCovariates), length(indicators));
	if (length(dataValues) - parametersEstimated < 5) {
		warning("Not enough degrees of freedom for accurate estimation of the basis coefficients")
		basis = splines::bs(dataValues, degree = 2, intercept = FALSE)
		degree = 2
		numKnots = 0
	} else if (length(dataValues) - parametersEstimated < 7) {
		basis = splines::bs(dataValues, degree = 2, knots = median(dataValues), intercept = FALSE)
		degree = 2
		numKnots = 1
	} else if (length(dataValues) - parametersEstimated < 10) {
		basis = splines::bs(dataValues, degree = 3, knots = median(dataValues), intercept = FALSE)
		degree = 3
		numKnots = 1
	} else if (length(dataValues) - parametersEstimated < 15) {
		basis = splines::bs(dataValues, degree = 3, knots = as.vector(quantile(dataValues)[2:4]), intercept = FALSE)
		degree = 3
		numKnots = 3
	} else {
		n = floor((length(dataValues) - parametersEstimated) / 6)
		range = max(dataValues) - min(dataValues)
		potentialKnotLocations = as.vector(quantile(dataValues, probs = seq(from = 1/n, to = 1 - 1/n, by = 1/n)))
		toRemove = integer(0)
		for (i in 2:length(potentialKnotLocations)) {
			j = i;
			while ((potentialKnotLocations[j] - potentialKnotLocations[i - 1]) < (range / 20)) {
				toRemove = c(toRemove, j) 
				j = j + 1
			}
		}
		if (length(toRemove) > 0) {
			potentialKnotLocations = potentialKnotLocations[-toRemove];
		}
		basis = splines::bs(dataValues, degree = 3, knots = potentialKnotLocations)
		degree = 3
		numKnots = length(potentialKnotLocations)
	}
	report = paste("Fit B-Spline of degree", degree, "with", numKnots, "internal knots to ", covariate, " data in the ", group, "group(s).")
	print(report)
	return(basis)
}
