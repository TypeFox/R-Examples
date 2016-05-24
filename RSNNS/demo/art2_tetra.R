library(RSNNS)

data(snnsData)
patterns <- snnsData$art2_tetra_low.pat

model <- art2(patterns, f2Units=5, learnFuncParams=c(0.99, 20, 20, 0.1, 0), updateFuncParams=c(0.99, 20, 20, 0.1, 0))

model

testPatterns <- snnsData$art2_tetra_high.pat
predictions <- predict(model, testPatterns)

library(scatterplot3d)

par(mfrow=c(2,2))
scatterplot3d(patterns, pch=encodeClassLabels(model$fitted.values))
scatterplot3d(testPatterns, pch=encodeClassLabels(predictions))
