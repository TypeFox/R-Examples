limmaPLM <- function(dataObject, intercept = TRUE, indicators = as.character(unique(dataObject$sampleInfo[,2])[-1]), continuousCovariates = NULL, groups = as.character(unique(dataObject$sampleInfo[,2])[-1]), groupFunctions = rep("AdditiveSpline", length(groups)), fitSplineFromData = TRUE, splineDegrees = rep(3, length(groups)), splineKnots = rep(0, length(groups)), splineKnotSpread = "quantile", ...) {
	
	if ((length(indicators) > 0) && (!(indicators %in% dataObject$sampleInfo[,2]))) {
		stop("invalid indicator variable")
	}
	
	## format the expression levels:
	expression = as.matrix(dataObject$expressionValues)
	rownames(expression) = dataObject$genes
	
	## stores DesignMatrix for call to lm():
	DesignMatrix = matrix(data = rep(1, length(dataObject$sampleInfo[, 1])), ncol = 1)
	if (length(indicators) > 0) { 
		for (condition in indicators) { 
			## create indicator variables:
			condition = toString(condition)
			indicator = as.numeric(dataObject$sampleInfo[,2] == condition)
			## incorporate indicator into design matrix:
			DesignMatrix = cbind(DesignMatrix, indicator)
			colnames(DesignMatrix)[length(DesignMatrix[1,])] = paste("Indicator.", condition, sep = "")
 		}
	}
	
	## Check for the inclusion of intercept term in model:
	if (!intercept) {
		DesignMatrix = DesignMatrix[ , -1]
	} else {
		colnames(DesignMatrix)[1] = "Intercept"
	}
	if ((!intercept) && (length(indicators) == 0)) {
		warning("No indicator variables or intercept.")
	}
	if (length(continuousCovariates > 0)) {
		for (covariate in continuousCovariates) {
			## get covariate:
			covariate = toString(covariate)
			covariateValues = dataObject$sampleInfo[,covariate]
			## coerce NA values to zero:
			covariateValues[is.na(covariateValues)] = 0
			if (length(groups) > 0) {
				## identify the different groups:					
				if (length(groups) != length(groupFunctions)) {
					stop("groups must have same length as groupFunctions")
				}
				GroupsSharingLinearFunction = character(0)
				GroupsSharingSplineFunction = character(0)
				for (i in 1:length(groups)) {
					group = groups[i]
					groupLocation = grep(group, dataObject$sampleInfo[,2])
					if (length(groupLocation) == 0) {
						stop("Invalid group name")
					}
					if (groupFunctions[i] == "AdditiveLinear") {
						additiveLinearCovariate = rep(0, length(dataObject$sampleInfo[,2]))
						additiveLinearCovariate[groupLocation] = covariateValues[groupLocation]
						DesignMatrix = cbind(DesignMatrix, additiveLinearCovariate)
						covarName = paste(covariate, ".", group, sep = "")
						colnames(DesignMatrix)[length(DesignMatrix[1, ])] = covarName
					} else if (groupFunctions[i] == "CommonLinear") {									GroupsSharingLinearFunction = c(GroupsSharingLinearFunction, group)
					} else if (groupFunctions[i] == "CommonSpline") {
						GroupsSharingSplineFunction = c(GroupsSharingSplineFunction, group)
					} else if (groupFunctions[i] == "AdditiveSpline") {
						if (sum(abs(covariateValues[groupLocation]) > 0 ) == 0) {
							## all covariate values in this group are 0
							stop("Cannot fit spline to group whose only covariate value is zero")
						}
						## fit B spline:
						if (fitSplineFromData) {
							splineBasis = fitBspline(covariateValues[groupLocation], continuousCovariates, indicators, group, covariate); 
						} else {
							numKnots = splineKnots[i];
							if (numKnots == 0) {
								KnotLocations = NULL
							} else {
								if (splineKnotSpread == "uniform") {
									KnotInterval = (max(covariateValues[groupLocation]) - min(covariateValues[groupLocation])) / (numKnots + 1)
								KnotLocations = seq(min(covariateValues[groupLocation]), max(covariateValues[groupLocation]), by = KnotInterval)[-c(1, numKnots + 2)]	
								} else if (splineKnotSpread == "quantile") {
									KnotLocations = as.vector(quantile(covariateValues[groupLocation], seq(0, 1, by = 1 / (numKnots + 1))))
									KnotLocations = KnotLocations[-c(1, numKnots + 2)]
								} else {
									stop("splineKnotSpread must be one of: uniform, quantile")
								}
							}
							splineBasis = splines::bs(covariateValues[groupLocation], knots = KnotLocations, degree = splineDegrees[i], intercept = FALSE)
						}
						OverallBasis = matrix(rep(0, length(covariateValues) * length(splineBasis[1,])), ncol = length(splineBasis[1,]))
						OverallBasis[groupLocation, ] = splineBasis
						for (k in 1:length(OverallBasis[1,])) {
							covarName = paste(group, "BasisFunction.", covariate, ".", k, sep = "")
							DesignMatrix = cbind(DesignMatrix, OverallBasis[, k])
							colnames(DesignMatrix)[length(DesignMatrix[1,])] = covarName
						} 	
					} else {
						stop("groupFunctions must be one of: CommonLinear, CommonSpline, AdditiveLinear, AdditiveSpline")
					}
				}
				## fit common spline to groups that share it:
				groupLocation = integer(0)
				commongroup = character(0)
				for (group in GroupsSharingSplineFunction) {
					groupLocation = c(groupLocation, grep(group, dataObject$sampleInfo[,2]))
					commongroup = c(commongroup, group)
				}
				if (length(groupLocation) > 0) {
					## fit B spline:
					if (sum(abs(covariateValues[groupLocation]) > 0 ) == 0) {
						## all covariate values in this group are 0
						stop("Cannot fit spline to group whose only covariate value is zero")
					}
					if (fitSplineFromData) {
						splineBasis = fitBspline(covariateValues[groupLocation], continuousCovariates, indicators, commongroup, covariate); 
					} else {
						numKnots = splineKnots[i];
						if (numKnots == 0) {
							KnotLocations = NULL
						} else {
							if (splineKnotSpread == "uniform") {
								KnotInterval = (max(covariateValues[groupLocation]) - min(covariateValues[groupLocation])) / (numKnots + 1)
								KnotLocations = seq(min(covariateValues[groupLocation]), max(covariateValues[groupLocation]), by = KnotInterval)[-c(1, numKnots + 2)]	
							} else if (splineKnotSpread == "quantile") {
								KnotLocations = as.vector(quantile(covariateValues[groupLocation], seq(0, 1, by = 1 / (numKnots + 1))))
								KnotLocations = KnotLocations[-c(1, numKnots + 2)]
							} else {
								stop("splineKnotSpread must be one of: uniform, quantile")
							}
						}
						splineBasis = splines::bs(covariateValues[groupLocation], knots = KnotLocations, degree = splineDegrees[i], intercept = FALSE)
					}					
					OverallBasis = matrix(rep(0, length(covariateValues) * length(splineBasis[1,])), ncol = length(splineBasis[1,]))
					OverallBasis[groupLocation, ] = splineBasis
					for (k in 1:length(OverallBasis[1,])) {
						covarName = paste("commonBasisFunction.", covariate, ".", k, sep = "")
						DesignMatrix = cbind(DesignMatrix, OverallBasis[ ,k])
						colnames(DesignMatrix)[length(DesignMatrix[1,])] = covarName
					}
				} 		
				## fit common linear function to groups that share it:
				groupLocation = integer(0)
				for (group in GroupsSharingLinearFunction) {
					groupLocation = c(groupLocation, grep(group, dataObject$sampleInfo[,2]))
				}
				if (length(groupLocation) > 0) {
					additiveLinearCovariate = rep(0, length(dataObject$sampleInfo[,2]))
					additiveLinearCovariate[groupLocation] = covariateValues[groupLocation]
					covarName = paste(covariate, ".Common", sep = "")
					DesignMatrix = cbind(DesignMatrix, additiveLinearCovariate)
						colnames(DesignMatrix)[length(DesignMatrix[1,])] = covarName
				}
			} else { # length of groups is zero
				stop("Each covariate requires specification of group and function fit")
			}
		}
	}
	## Check if there is a Design Matrix:
	if ((dim(DesignMatrix) == 0) || (is.null(DesignMatrix))) {
		stop("model contains no design matrix")
	}
    MArrayLMobject = limma::lmFit(expression, design = DesignMatrix, ...)
	return(MArrayLMobject)
}
