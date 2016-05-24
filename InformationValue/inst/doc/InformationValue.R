## ---- echo=FALSE, results='asis', fig.show='hold', message=FALSE, warning=FALSE----
library(InformationValue)

## ---- results='asis', fig.show='hold', fig.height=5, fig.width=8, prompt=TRUE, highlight=TRUE, collapse=TRUE, message=FALSE, warning=FALSE----
data('ActualsAndScores')
plotROC(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)

## ---- echo=TRUE, results='hide', fig.show='hide', prompt=TRUE, highlight=TRUE, tidy=TRUE, highlight=TRUE, collapse=TRUE----
sensMat <- plotROC(actuals=ActualsAndScores$Actuals,  predictedScores=ActualsAndScores$PredictedScores, returnSensitivityMat = TRUE)

## ---- results='asis', prompt=TRUE, highlight=TRUE, collapse=TRUE, message=FALSE, warning=FALSE, highlight=TRUE, tidy=TRUE----
sensitivity(actuals = ActualsAndScores$Actuals, predictedScores = ActualsAndScores$PredictedScores)

## ---- results='asis', prompt=TRUE, highlight=TRUE, message=FALSE, warning=FALSE, highlight=TRUE, tidy=TRUE----
max_sens_cutoff <- optimalCutoff(actuals=ActualsAndScores$Actuals, predictedScores = ActualsAndScores$PredictedScores, optimiseFor='Ones')  # determine cutoff to maximize sensitivity.

print(max_sens_cutoff)  # This would be cut-off score that achieved maximum sensitivity.

sensitivity(actuals = ActualsAndScores$Actuals, predictedScores = ActualsAndScores$PredictedScores, threshold=max_sens_cutoff)

## ---- results='asis', prompt=TRUE, highlight=TRUE, collapse=FALSE, message=FALSE, warning=FALSE, highlight=TRUE, tidy=TRUE----
specificity(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)

## ---- results='asis', prompt=TRUE, highlight=TRUE, collapse=FALSE, message=FALSE, warning=FALSE, tidy=TRUE----
specificity(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores, threshold = 0.35)

## ---- results='asis', prompt=TRUE, highlight=TRUE, collapse=FALSE, message=FALSE, warning=FALSE, highlight=TRUE, tidy=TRUE----
precision(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)

## ---- results='asis', prompt=TRUE, highlight=TRUE, collapse=FALSE, message=FALSE, warning=FALSE, highlight=TRUE, tidy=TRUE----
npv(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)

## ---- results='asis', prompt=TRUE, highlight=TRUE, collapse=TRUE, message=FALSE, warning=FALSE, tidy=TRUE, comment='#'----
youdensIndex(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)

## ---- results='asis', prompt=TRUE, highlight=TRUE, collapse=TRUE, message=FALSE, warning=FALSE, tidy=TRUE, comment='#'----
misClassError(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores, threshold=0.5)

## ---- results='asis', prompt=TRUE, highlight=TRUE, collapse=TRUE, message=FALSE, warning=FALSE, tidy=TRUE, comment='#'----
kappaCohen(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)

## ---- results='asis', prompt=TRUE, highlight=TRUE, collapse=TRUE, message=FALSE, warning=FALSE, tidy=TRUE, comment='#'----
Concordance(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)

## ---- results='asis', prompt=TRUE, highlight=TRUE, collapse=TRUE, message=FALSE, warning=FALSE, tidy=TRUE, comment='#'----
somersD(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)

## ---- results='asis', prompt=TRUE, highlight=TRUE, collapse=TRUE, message=FALSE, warning=FALSE, tidy=TRUE, comment='#'----
ks_stat(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)

## ---- , results='asis', fig.show='hold', fig.height=5, fig.width=8, prompt=TRUE, highlight=TRUE, collapse=TRUE, message=FALSE, warning=FALSE----
ks_plot(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)

## ---- results='asis', prompt=TRUE, highlight=TRUE, collapse=TRUE, message=FALSE, warning=FALSE, tidy=TRUE, comment='#'----
optimalCutoff(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores)  # returns cutoff that gives minimum misclassification error.

## ---- results='asis', prompt=TRUE, highlight=TRUE, collapse=TRUE, message=FALSE, warning=FALSE, tidy=TRUE, comment='#'----
optimalCutoff(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores, optimiseFor="Both")  # returns cutoff that gives maximum of Youden's J Index

## ---- results='asis', prompt=TRUE, highlight=TRUE, collapse=TRUE, message=FALSE, warning=FALSE, tidy=TRUE, comment='#'----
sens_table <- optimalCutoff(actuals=ActualsAndScores$Actuals, predictedScores=ActualsAndScores$PredictedScores, optimiseFor="Both", returnDiagnostics=TRUE)$sensitivityTable

## ---- results='hide', prompt=TRUE, highlight=TRUE, collapse=TRUE, message=FALSE, warning=FALSE, tidy=TRUE, comment='#'----
WOE(X=SimData$X.Cat, Y=SimData$Y.Binary)

## ---- results='hide', prompt=TRUE, highlight=TRUE, collapse=TRUE, message=FALSE, warning=FALSE, tidy=TRUE, comment='#'----
options(scipen = 999, digits = 2)
WOETable(X=SimData$X.Cat, Y=SimData$Y.Binary)

## ---- results='asis', prompt=TRUE, highlight=TRUE, collapse=TRUE, message=FALSE, warning=FALSE, tidy=TRUE, comment='#'----
options(scipen = 999, digits = 4)
IV(X=SimData$X.Cat, Y=SimData$Y.Binary)

