plot.DEresults <- function(x, covariate, geneNumber = 1, plmDEobject, loess = TRUE, legend = TRUE, legend.coor = "topright",  ... ) {
	Indicator.full = x$modelForm.fullModel[[1]]
	quantitativeCovariateList = x$modelForm.fullModel[[2]]
	fullform = as.formula(paste("unlist(plmDEobject$expressionValues[geneNumber, ]) ~ ", x$PredictorFormula.fullModel))	
	
	## fit model to selected gene:
	if (!(x$GLMinfo[[1]][1])) {
		fullModel = lm(fullform)
	} else if (x$GLMinfo[["NegativeBinomialUnknownDispersion"]][1]) {
		fullModel = MASS::glm.nb(fullform)
	} else {
		fullModel = glm(fullform, family = x$GLMinfo[["family"]][1])
	}
	fittedValues.full = fullModel$fitted.values
	
	## using indicator, assign 0 to samples with missing 
	## measurements of covariate:
	covariateMeasurements = plmDEobject$sampleInfo[, covariate]
	measurementsMissing = is.na(covariateMeasurements)
	covariateMeasurements[measurementsMissing] = 0 
	## Plot the fit of the model:
	plot(covariateMeasurements, plmDEobject$expressionValues[geneNumber, ], type = "n", xlab = covariate, ylab = paste("Expression Level of gene",geneNumber,"across samples"), ...)
	nonmissingMeasurements = covariateMeasurements[!(measurementsMissing | covariateMeasurements == 0)]
	connectable = order(nonmissingMeasurements);
	CorrespondingfittedValues.full = fittedValues.full[!(measurementsMissing | covariateMeasurements == 0)]
	connectable = order(nonmissingMeasurements);
	lines(nonmissingMeasurements[connectable], CorrespondingfittedValues.full[connectable], col = 2, lwd = 1.5)
	points(nonmissingMeasurements, CorrespondingfittedValues.full, col = 2, pch = 16, cex = 0.75)	
	xgrid = seq(min(covariateMeasurements), max(covariateMeasurements), length.out = 20)
	## loess fit of expression level based on covariate measure:
	lo = loess(unlist(plmDEobject$expressionValues[geneNumber, ]) ~ plmDEobject$sampleInfo[, covariate])
	lines(xgrid, predict(lo, xgrid), lty = 3)
	
	## add points to plot based on their group:
	groups = unique(plmDEobject$sampleInfo[,2]);
	plottingChar = 1;
	for (group in groups) {
		group = as.character(group)
		location = grep(group, as.character(plmDEobject$sampleInfo[,2]))
		points(covariateMeasurements[location], unlist(plmDEobject$expressionValues[geneNumber, ])[location], pch = plottingChar)
		plottingChar = plottingChar + 1;		
	}
	
	## check for reduced model:
	if (length(x$PredictorFormula.reducedModel) > 0) {
		Indicator.reduced = x$modelForm.reducedModel[[1]]
		quantitativeCovariateList.reduced = x$modelForm.reducedModel[[2]]
		reducedform = as.formula(paste("unlist(plmDEobject$expressionValues[geneNumber, ]) ~ ", x$PredictorFormula.reducedModel))
		## fit reduced model to selected gene:
		if (!(x$GLMinfo[[1]][1])) {
			reducedModel = lm(reducedform)
		} else if (x$GLMinfo[["NegativeBinomialUnknownDispersion"]][1]) {
			reducedModel = MASS::glm.nb(reducedform)
		} else {
			reducedModel = glm(reducedform, family = x$GLMinfo[["family"]][1])
		}
		fittedVals.reduced = reducedModel$fitted.values
		CorrespondingfittedValues.reduced = fittedVals.reduced[!(measurementsMissing | covariateMeasurements == 0)]
		lines(nonmissingMeasurements[connectable], CorrespondingfittedValues.reduced[connectable], col = 4, lwd = 1, lty = 2)
		points(nonmissingMeasurements, CorrespondingfittedValues.reduced, col = 4, pch = 16, cex = 0.55)
    }
    if (legend == TRUE) {
        legend(legend.coor, legend = groups, pch = c(1:length(groups)))
    }
}
