library(RSNNS)

data(snnsData)

patterns <- snnsData$art1_letters.pat

model <- assoz(patterns, dimX=7, dimY=5)

#model$fitted.values

actMaps <- matrixToActMapList(model$fitted.values, nrow=7)

par(mfrow=c(3,3))
for (i in 1:9) plotActMap(actMaps[[i]])

#predict(model, patterns)