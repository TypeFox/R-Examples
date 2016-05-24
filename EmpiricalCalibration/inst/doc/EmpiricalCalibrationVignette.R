## ---- echo = FALSE, message = FALSE, warning = FALSE---------------------
library(EmpiricalCalibration)

## ------------------------------------------------------------------------
data(sccs)
drugOfInterest <- sccs[sccs$groundTruth == 1, ]
drugOfInterest
exp(drugOfInterest$logRr)
computeTraditionalP(drugOfInterest$logRr, drugOfInterest$seLogRr)

## ------------------------------------------------------------------------
data(sccs)
negatives <- sccs[sccs$groundTruth == 0, ]
head(negatives)

## ------------------------------------------------------------------------
plotForest(negatives$logRr, negatives$seLogRr, negatives$drugName)

## ------------------------------------------------------------------------
null <- fitNull(negatives$logRr, negatives$seLogRr)
null

## ------------------------------------------------------------------------
plotCalibration(negatives$logRr,negatives$seLogRr)

## ------------------------------------------------------------------------
plotCalibrationEffect(negatives$logRr,negatives$seLogRr, null = null)

## ------------------------------------------------------------------------
p <- calibrateP(null, drugOfInterest$logRr, drugOfInterest$seLogRr)
p

## ------------------------------------------------------------------------
plotCalibrationEffect(negatives$logRr,
                      negatives$seLogRr, 
                      drugOfInterest$logRr, 
                      drugOfInterest$seLogRr, 
                      null)

## ------------------------------------------------------------------------
null <- fitMcmcNull(negatives$logRr, negatives$seLogRr)
null

## ------------------------------------------------------------------------
plotMcmcTrace(null)

## ------------------------------------------------------------------------
p <- calibrateP(null, drugOfInterest$logRr, drugOfInterest$seLogRr)
p

## ----tidy=TRUE,evale=TRUE------------------------------------------------
citation("EmpiricalCalibration")

