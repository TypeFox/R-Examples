fitGAPLM <- function(dataObject, generalizedLM = FALSE, family = poisson(link = log), NegativeBinomialUnknownDispersion = FALSE, test = "LRT", weights = NULL, offset = NULL, pValueAdjustment = "fdr", significanceLevel = 0.05, indicators.fullModel = as.character(unique(dataObject$sampleInfo[,2])[-1]), continuousCovariates.fullModel = NULL, groups.fullModel = as.character(unique(dataObject$sampleInfo[,2])[-1]), groupFunction.fullModel = rep("AdditiveSpline", length(groups.fullModel)), fitSplineFromData.fullModel = TRUE, splineDegrees.fullModel = rep(3, length(groups.fullModel)), splineKnots.fullModel = rep(0, length(groups.reducedModel)), compareToReducedModel = FALSE, indicators.reducedModel = as.character(unique(dataObject$sampleInfo[,2])[-1]), continuousCovariates.reducedModel = NULL, groups.reducedModel = as.character(unique(dataObject$sampleInfo[,2])[-1]), groupFunction.reducedModel = rep("AdditiveSpline", length(groups.reducedModel)), fitSplineFromData.reducedModel = TRUE, splineDegrees.reducedModel = rep(3, length(groups.reducedModel)), splineKnots.reducedModel = rep(0, length(groups.reducedModel)), splineKnotSpread = "quantile") { 	
	
	if ((length(indicators.fullModel) > 0) && (!(indicators.fullModel %in% dataObject$sampleInfo[,2]))) {
		stop("invalid indicator variable")
	}
	if ((length(indicators.reducedModel) > 0) && (!(indicators.reducedModel %in% dataObject$sampleInfo[,2]))) {
		stop("invalid indicator variable")
	}
	## format the expression levels:
	expression = as.matrix(dataObject$expressionValues)
	## stores formula for call to lm();
	stringFormula.full = ""
    Indicator.full = list()
	if (length(indicators.fullModel) > 0) { 
		for (condition in indicators.fullModel) { 
			## create indicator variables:
			condition = toString(condition)
			Indicator.full[[condition]] = as.numeric(dataObject$sampleInfo[,2] == condition)
			## incorporate indicator into formula:
			stringFormula.full = paste(stringFormula.full, "  Indicator.full[[\'", condition, "\']] + ", sep = "")
		}
	}
	quantitativeCovariateList = list()
	if (length(continuousCovariates.fullModel > 0)) {
		for (covariate in continuousCovariates.fullModel) {
			## get covariate:
			covariate = toString(covariate)
			covariateValues = dataObject$sampleInfo[,covariate]
			## coerce NA values to zero:
			covariateValues[is.na(covariateValues)] = 0
			if (length(groups.fullModel) > 0) {
				## identify the different groups:					
				if (length(groups.fullModel) != length(groupFunction.fullModel)) {
					stop("groups.fullModel must have same length as groupFunction.fullModel")
				}
				GroupsSharingLinearFunction = character(0)
				GroupsSharingSplineFunction = character(0)
				for (i in 1:length(groups.fullModel)) {
					group = groups.fullModel[i]
					groupLocation = grep(group, dataObject$sampleInfo[,2])
					if (length(groupLocation) == 0) {
						stop("Invalid group name")
					}
					if (groupFunction.fullModel[i] == "AdditiveLinear") {
						additiveLinearCovariate = rep(0, length(dataObject$sampleInfo[,2]))
						additiveLinearCovariate[groupLocation] = covariateValues[groupLocation]
						covarName = paste(covariate, ".", group, sep = "")
						quantitativeCovariateList[[covarName]] = additiveLinearCovariate
					} else if (groupFunction.fullModel[i] == "CommonLinear") {							GroupsSharingLinearFunction = c(GroupsSharingLinearFunction, group)
					} else if (groupFunction.fullModel[i] == "CommonSpline") {
						GroupsSharingSplineFunction = c(GroupsSharingSplineFunction, group)
					} else if (groupFunction.fullModel[i] == "AdditiveSpline") {
						if (sum(abs(covariateValues[groupLocation]) > 0 ) == 0) {
							## all covariate values in this group are 0
							stop("Cannot fit spline to group whose only covariate value is zero")
						}
						## fit B spline:
						if (fitSplineFromData.fullModel) {
							splineBasis = fitBspline(covariateValues[groupLocation], continuousCovariates.fullModel, indicators.fullModel, group, covariate); 
						} else {
							numKnots = splineKnots.fullModel[i];
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
							splineBasis = splines::bs(covariateValues[groupLocation], knots = KnotLocations, degree = splineDegrees.fullModel[i], intercept = FALSE)
						}
						OverallBasis = matrix(rep(0, length(covariateValues) * length(splineBasis[1,])), ncol = length(splineBasis[1,]), nrow = length(covariateValues))
						OverallBasis[groupLocation, ] = splineBasis
						for (k in 1:length(OverallBasis[1,])) {
							covarName = paste(group, "BasisFunction", k, ".", covariate, sep = "")
							quantitativeCovariateList[[covarName]] = OverallBasis[,k]
						} 	
					} else {
						stop("groupFunction must be one of: CommonLinear, CommonSpline, AdditiveLinear, AdditiveSpline")
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
					if (fitSplineFromData.fullModel) {
						splineBasis = fitBspline(covariateValues[groupLocation], continuousCovariates.fullModel, indicators.fullModel, commongroup, covariate); 
					} else {
						numKnots = splineKnots.fullModel[i];
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
						splineBasis = splines::bs(covariateValues[groupLocation], knots = KnotLocations, degree = splineDegrees.fullModel[i], intercept = FALSE)
					}					
					OverallBasis = matrix(rep(0, length(covariateValues) * length(splineBasis[1,])), ncol = length(splineBasis[1,]), nrow = length(covariateValues))
					OverallBasis[groupLocation, ] = splineBasis
					for (k in 1:length(OverallBasis[1,])) {
						covarName = paste("commonBasisFunction", k, ".", covariate, sep = "")
						quantitativeCovariateList[[covarName]] = OverallBasis[,k]
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
					covarName = paste(covariate, ".common", sep = "")
					quantitativeCovariateList[[covarName]] = additiveLinearCovariate
				}
			} else { # length of groups.fullModel is zero
				stop("Each covariate requires specification of group and function fit")
			}
		}
		for (covarName in names(quantitativeCovariateList)) { 
			## add continuous covariates to formula for lm:
			stringFormula.full = paste(stringFormula.full, " quantitativeCovariateList[[\'", covarName, "\']] + ", sep = "")
		}
	}
	## Check if there is a model to fit:
	if (nchar(stringFormula.full) == 0) {
		stop("full model contains no indicators or numerical covariates")
	}
	stringFormula.full = R.oo::trim(stringFormula.full)
	length = nchar(stringFormula.full)
	if (substring(stringFormula.full, length, length) == "+") {
		stringFormula.full = substring(stringFormula.full, 1, length - 1)
	}	
	## Check if Reduced Model is requested:
	if (!compareToReducedModel) { # test all coefficients for significance
		PvaluesForCoefficientsSignificance = rep(0, length(dataObject$genes))
		if (generalizedLM) { # fit generalized linear model
			if (NegativeBinomialUnknownDispersion) {
				SignificanceReport = "Pr(Chi)"
				for (j in 1:length(dataObject$genes)) {
				PvaluesForCoefficientsSignificance[j] = anova(MASS::glm.nb(as.formula(paste("expression[j, ] ~ ", stringFormula.full)), weights = weights), MASS::glm.nb(expression[j, ] ~ 1, weights = weights), test = test)[[SignificanceReport]][2]
				}		
			} else {
				whichTest = grep(test, c("ChiSq", "LRT", "Rao", "F"))
				if (whichTest < 4) {
					SignificanceReport = "Pr(>Chi)"
				} else if (whichTest == 4) {
					SignificanceReport = "Pr(>F)"
				} else {
					stop("test must be one of: ChiSq, LRT, Rao, F")
				}
				for (j in 1:length(dataObject$genes)) {
					PvaluesForCoefficientsSignificance[j] = anova(glm(as.formula(paste("expression[j, ] ~ ", stringFormula.full)), family = family, weights = weights, offset = offset), glm(expression[j, ] ~ 1, weights = weights, offset = offset), test = test)[[SignificanceReport]][2]
				}
			}
		} else { # fit linear model and test for significance
			for (j in 1:length(dataObject$genes)) {
				PvaluesForCoefficientsSignificance[j] = 1 - pf(summary(lm(as.formula(paste("expression[j, ] ~ ", stringFormula.full)), weights, offset))$fstatistic[1], summary(lm(as.formula(paste("expression[j, ] ~ ", stringFormula.full)), weights = weights, offset = offset))$fstatistic[2], summary(lm(as.formula(paste("expression[j, ] ~ ", stringFormula.full)), weights = weights, offset = offset))$fstatistic[3])		
			}
		}
		# data frame of genes and p-values
		DifferentiallyExpressed = data.frame(gene = dataObject$genes, p.val = PvaluesForCoefficientsSignificance)
		stringFormula.reduced = NULL
		
	} else { # Construct another formula for Reduced Model:
        Indicator.reduced = list()
		stringFormula.reduced = ""
		if (length(indicators.reducedModel > 0)) { 
			for (condition in indicators.reducedModel) { 
				## create indicator variables:
				condition = toString(condition)
				Indicator.reduced[[condition]] = as.numeric(dataObject$sampleInfo[,2] == condition)
				## incorporate indicator into formula:
				stringFormula.reduced = paste(stringFormula.reduced, " Indicator.reduced[[\'", condition, "\']] + ", sep = "")
			}
		}
		quantitativeCovariateList.reduced = list()
		if (length(continuousCovariates.reducedModel > 0)) {
			for (covariate in continuousCovariates.reducedModel) {
				## get covariate:
				covariate = toString(covariate)
				covariateValues = dataObject$sampleInfo[,covariate]
				## coerce NA values to zero:
				covariateValues[is.na(covariateValues)] = 0
				if (length(groups.reducedModel) > 0) {
					## identify the different groups:					
					if (length(groups.reducedModel) != length(groupFunction.reducedModel)) {
						stop("groups.reducedModel must have same length as groupFunction.reducedModel")
					}
					GroupsSharingLinearFunction = character(0)
					GroupsSharingSplineFunction = character(0)
					for (i in 1:length(groups.reducedModel)) {
						group = groups.reducedModel[i]
						groupLocation = grep(group, dataObject$sampleInfo[,2])
						if (length(groupLocation) == 0) {
							stop("Invalid group name")
						}
						if (groupFunction.reducedModel[i] == "AdditiveLinear") {
							additiveLinearCovariate = rep(0, length(dataObject$sampleInfo[,2]))
							additiveLinearCovariate[groupLocation] = covariateValues[groupLocation]
							covarName = paste(covariate, ".", group, sep = "")
							quantitativeCovariateList.reduced[[covarName]] = additiveLinearCovariate
						} else if (groupFunction.reducedModel[i] == "CommonLinear") {		
							GroupsSharingLinearFunction = c(GroupsSharingLinearFunction, group)
						} else if (groupFunction.reducedModel[i] == "CommonSpline") {
							GroupsSharingSplineFunction = c(GroupsSharingSplineFunction, group)
						} else if (groupFunction.reducedModel[i] == "AdditiveSpline") {
							if (sum(abs(covariateValues[groupLocation]) > 0 ) == 0) {
								## all covariate values in this group are 0
								stop("Cannot fit spline to group whose only covariate value is zero")
							}
							## fit B spline:
						if (fitSplineFromData.reducedModel) {
							splineBasis = fitBspline(covariateValues[groupLocation], continuousCovariates.reducedModel, indicators.reducedModel, group, covariate); 
						} else {
							numKnots = splineKnots.reducedModel[i];
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
							splineBasis = splines::bs(covariateValues[groupLocation], knots = KnotLocations, degree = splineDegrees.reducedModel[i], intercept = FALSE)
						}
							OverallBasis = matrix(rep(0, length(covariateValues) * length(splineBasis[1,])), ncol = length(splineBasis[1,]), nrow = length(covariateValues))
							OverallBasis[groupLocation, ] = splineBasis
							for (k in 1:length(OverallBasis[1,])) {
								covarName = paste(group, "BasisFunction", k, ".", covariate, sep = "")
								quantitativeCovariateList.reduced[[covarName]] = OverallBasis[,k]
							} 	
						} else {
							stop("groupFunction must be one of: CommonLinear, CommonSpline, AdditiveLinear, AdditiveSpline")
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
							## fit B spline:
						if (fitSplineFromData.reducedModel) {
							splineBasis = fitBspline(covariateValues[groupLocation], continuousCovariates.reducedModel, indicators.reducedModel, group, covariate);
						} else {
							numKnots = splineKnots.reducedModel[i];
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
							splineBasis = splines::bs(covariateValues[groupLocation], knots = KnotLocations, degree = splineDegrees.reducedModel[i], intercept = FALSE)
						}
						OverallBasis = matrix(rep(0, length(covariateValues) * length(splineBasis[1,])), ncol = length(splineBasis[1,]), nrow = length(covariateValues))
						OverallBasis[groupLocation, ] = splineBasis
						for (k in 1:length(OverallBasis[1,])) {
							covarName = paste("commonBasisFunction", k, ".", covariate, sep = "")
							quantitativeCovariateList.reduced[[covarName]] = OverallBasis[,k]
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
						covarName = paste(covariate, ".common", sep = "")
						quantitativeCovariateList.reduced[[covarName]] = additiveLinearCovariate
					}
				} else { # length of groups.reducedModel is zero
					stop("Each covariate requires specification of group and function fit")
				}
			}
			for (covarName in names(quantitativeCovariateList.reduced)) { 
				## add continuous covariates to formula for lm:
				stringFormula.reduced = paste(stringFormula.reduced, " quantitativeCovariateList.reduced[[\'", covarName, "\']] + ", sep = "")
			}
		}
		## Ensure model is not empty:
		if (nchar(stringFormula.reduced) == 0) {
			stop("reduced model has no indicators or numerical covariates")		
		}
		stringFormula.reduced = R.oo::trim(stringFormula.reduced)
		length = nchar(stringFormula.reduced)
		if (substring(stringFormula.reduced, length, length) == "+") {
			stringFormula.reduced = substring(stringFormula.reduced, 1, length - 1)
		}	
				
		## Comparing the full and reduced models
		PvaluesForSignificanceOfFullModel = rep(0, length(dataObject$genes))
		if (generalizedLM) { # fit generalized linear model
			if (NegativeBinomialUnknownDispersion) {
				SignificanceReport = "Pr(Chi)"
				for (j in 1:length(dataObject$genes)) {
				PvaluesForSignificanceOfFullModel[j] = anova(MASS::glm.nb(as.formula(paste("expression[j, ] ~ ", stringFormula.full)), weights = weights), MASS::glm.nb(as.formula(paste("expression[j, ] ~ ", stringFormula.reduced)), weights = weights), test = test)[[SignificanceReport]][2]
				}
			} else {
				whichTest = grep(test, c("ChiSq", "LRT", "Rao", "F"))
				if (whichTest < 4) {
					SignificanceReport = "Pr(>Chi)"
				} else if (whichTest == 4) {
					SignificanceReport = "Pr(>F)"
				} else {
					stop("test must be one of: ChiSq, LRT, Rao, F")
				}
				for (j in 1:length(dataObject$genes)) {
					PvaluesForSignificanceOfFullModel[j] = anova(glm(as.formula(paste("expression[j, ] ~ ", stringFormula.full)), family = family, weights = weights, offset = offset), glm(as.formula(paste("expression[j, ] ~ ", stringFormula.reduced)), weights = weights, offset = offset), test = test)[[SignificanceReport]][2]
				}
			}
		} else { # use linear model
			for (j in 1:length(dataObject$genes)) {	
				PvaluesForSignificanceOfFullModel[j] = anova(lm(as.formula(paste("expression[j,] ~ ", stringFormula.full)), weights = weights, offset = offset), lm(as.formula(paste("expression[j,] ~ ", stringFormula.reduced)), weights = weights, offset = offset))$"Pr(>F)"[2]	
			}
		}
		DifferentiallyExpressed = data.frame(gene = dataObject$genes, p.val = PvaluesForSignificanceOfFullModel)
	}
	
	## Apply p-value cut-off to determine significance
	## of differential expression level and find DE genes.
	AdjustedP = p.adjust(DifferentiallyExpressed$p.val, method = pValueAdjustment)
	DifferentiallyExpressed$p.val.Adjusted = AdjustedP
	DifferentiallyExpressed$Significant = (AdjustedP < significanceLevel)
	DEgenes = data.frame(gene = DifferentiallyExpressed$gene[DifferentiallyExpressed$Significant], p.val.Adjusted = AdjustedP[DifferentiallyExpressed$Significant])
	
	## order DEgenes in terms of increasing p-value
	## so that most differentially expressed genes appear at the top
	DEgenes = DEgenes[order(DEgenes$p.val.Adjusted), ]
	row.names(DEgenes) = NULL;
	
	## keep track of the model's form:
	modelForm.fullModel = list(Indicator.full, quantitativeCovariateList)
	if (length(stringFormula.reduced) > 0) {
		modelForm.reducedModel = list(Indicator.reduced, quantitativeCovariateList.reduced)
	} else {
		modelForm.reducedModel = NULL
	}
	## keep track of whether model is generalized linear model:
	GLMinfo = list(generalizedLM, family, NegativeBinomialUnknownDispersion)
	
	## store results in class DEresults:
	results = list(allgenes = DifferentiallyExpressed, DEgenes = DEgenes, PredictorFormula.fullModel = stringFormula.full, PredictorFormula.reducedModel = stringFormula.reduced, modelForm.fullModel = modelForm.fullModel, modelForm.reducedModel = modelForm.reducedModel, GLMinfo = GLMinfo)
	class(results) = "DEresults"
	return(results)
}