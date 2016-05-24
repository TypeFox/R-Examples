## ----exampledata---------------------------------------------------------
library('vtreat')
dTrainC <- data.frame(x=c('a','a','a','b','b',NA),
                      y=c(FALSE,FALSE,TRUE,FALSE,TRUE,TRUE))
treatmentsC <- designTreatmentsC(dTrainC,colnames(dTrainC),'y',TRUE,
                                 verbose=FALSE)
dTrainCTreatedUnscaled <- prepare(treatmentsC,dTrainC,pruneSig=c(),scale=FALSE)
dTrainCTreatedScaled <- prepare(treatmentsC,dTrainC,pruneSig=c(),scale=TRUE)

## ----printorig-----------------------------------------------------------
print(dTrainC)

## ----printunscaled-------------------------------------------------------
print(dTrainCTreatedUnscaled)

## ----printscaled---------------------------------------------------------
print(dTrainCTreatedScaled)

## ----check---------------------------------------------------------------
slopeFrame <- data.frame(varName = treatmentsC$scoreFrame$varName,
                         stringsAsFactors = FALSE)
slopeFrame$mean <-
  vapply(dTrainCTreatedScaled[, slopeFrame$varName, drop = FALSE], mean,
         numeric(1))
slopeFrame$slope <- vapply(slopeFrame$varName,
                           function(c) {
                             lm(paste('y', c, sep = '~'),
                                data = dTrainCTreatedScaled)$coefficients[[2]]
                           },
                           numeric(1))
slopeFrame$sig <- vapply(slopeFrame$varName,
                         function(c) {
                           treatmentsC$scoreFrame[treatmentsC$scoreFrame$varName == c, 'sig']
                         },
                         numeric(1))
slopeFrame$badSlope <-
  ifelse(is.na(slopeFrame$slope), TRUE, abs(slopeFrame$slope - 1) > 1.e-8)
print(slopeFrame)

