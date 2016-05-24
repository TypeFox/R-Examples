library(RSNNS)

data(snnsData)

patterns <- snnsData$art1_letters.pat

inputMaps <- matrixToActMapList(patterns, nrow=7)
par(mfrow=c(3,3))
for (i in 1:9) plotActMap(inputMaps[[i]])

model <- art1(patterns, dimX=7, dimY=5)
#updateFunc="ART1_Synchronous"

encodeClassLabels(model$fitted.values)

summary(model)
