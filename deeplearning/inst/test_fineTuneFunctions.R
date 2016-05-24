# New a DArch instance
source("inst/dropout.R")


darch = newDArch(c(10,20,1), batchSize = 10)
setDropoutOneMaskPerEpoch(darch) = F
setFineTuneFunction(darch) <- minimizeClassifier
setFineTuneFunction(darch) <- backpropagation
setFineTuneFunction(darch) <- backpropSGD
darch = generateDropoutMasksForDarch(darch)

# New a dataset
input <- matrix(runif(250), 50, 5)
target <- rowSums(cos(input) + sin(input)^2)

mean_v <- mean(target)
target <- as.numeric(target > mean_v * 1.02 )

input_test <- matrix(runif(100), 20, 5)
target_test <- rowSums(cos(input_test) + sin(input_test)^2)
mean_v <- mean(target_test)
target_test <- as.numeric(target_test > mean_v * 1.02)

# Compare with the benchmark - backpropagation

darch_1 = darch( x = input,
               y = target,
               layers = c(5, 100, 50, 1),
               # darch = darch,
               darch.layerFunctionDefault = rectified_linear_unit_function,
               darch.layerFunctions = c("3" = sigmoidUnitDerivative),
               darch.bootstrap = F,
               darch.isBin = F,
               darch.isClass = F,
               darch.learnRateWeights = 0.01,
               darch.learnRateBiases = 0.01,
               darch.dropoutInput = 0.,
               darch.dropoutHidden = 0.,
               darch.fineTuneFunction = backpropagation, # finetune_SGD,
               darch.batchSize = 10,
               darch.numEpochs = 50
               )


darch_2 = darch( x = input,
               y = target,
               layers = c(5, 100, 50, 1),
               # darch = darch,
               darch.layerFunctionDefault = rectified_linear_unit_function,
               darch.layerFunctions = c("3" = sigmoidUnitDerivative),
               darch.bootstrap = F,
               darch.isBin = F,
               darch.isClass = F,
               darch.learnRateWeights = 0.01,
               darch.learnRateBiases = 0.01,
               darch.dropoutInput = 0.,
               darch.dropoutHidden = 0.,
               # darch.errorFunction = crossEntropyError,
               darch.fineTuneFunction = finetune_SGD_bn,
               errorFunc = meanSquareErr,
               darch.batchSize = 10,
               darch.numEpochs = 50
)

AR(darch_1)
AR(darch_2)

AR(darch_1, input_test, target_test)
AR(darch_2, input_test, target_test)

plot(predict(darch_1), predict(darch_2))

# Just use the finetuneDArch method.
# This function should be seperated to a train_dnn function
dataset <- createDataSet(input, target)

darch3 = fineTuneDArch(darch_1, dataset,
                       dataSetValid = NULL,
                       numEpochs = 5,
                       bootstrap = F,
                       isBin = T,
                       isClass = T,
                       stopErr = -Inf,
                       stopClassErr = -Inf,
                       stopValidErr = -Inf,
                       stopValidClassErr = 101
                       )

# Use the fineTune function directly

darch2 = darch

# Backpropagation/ Steepest Descent
darch2 = backpropagation(darch, dataset@data, dataset@targets)

# Conjugate Gradient Descent - Doesn't seem to work well. A bug in the code?
# darch3 = backpropCGD(darch, dataset@data, dataset@targets, length = 3, switchLayers = 0)

# Modified Steepest Gradient Descent
darch2 <- backpropSGD(darch, dataset@data, dataset@targets, crossEntropyErr)

testFunc2(darch2, dataset@data, dataset@targets, "Train set")

getLayer(darch,1)[[1]][1,]
getLayer(darch2,1)[[1]][1,]
getLayer(darch3,1)[[1]][1,]

testFunc2(darch3, dataset@data, dataset@targets, "Train Set")


gr1 <- calcGradient(par, darch2, dims, data, target, crossEntropyErr, epochSwitch)
gr2 <- fr(par, darch2, dims, data, target, epochSwitch )
gr1 - gr2

