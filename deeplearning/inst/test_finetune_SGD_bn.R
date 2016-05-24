# a toy model to test the finetune_SGD_bn function

input <- matrix(runif(100), 50, 2)
target <- rowSums(input + input^2)


# new a darch instance using new_darch
darch <- new_dnn(c(2, 20, 30, 20, 1))


for (i in 1:100) {
  darch <- generateDropoutMasksForDarch(darch, dropout_input = 0.2, dropout_hidden = 0.5)
  darch <- finetune_SGD_bn(darch, input, target,
                           learn_rate_weight = exp(-10),
                           learn_rate_bias = exp(-10),
                           learn_rate_gamma = exp(-10),
                           with_BN = T)
  darch <- calcualte_population_mu_sigma(darch, input)
  ret <- mseError(target, predict(darch, newdata = input))
  print(paste0(ret[[1]], ", ", ret[[2]]))
}


















darch = darch( x = input,
                 y = target,
                 layers = c(2, 10, 1),
                 # darch = darch,
                 darch.layerFunctionDefault = rectified_linear_unit_function,
                 darch.layerFunctions = c("2" = linearUnitDerivative),
                 darch.bootstrap = F,
                 darch.isBin = F,
                 darch.isClass = F,
                 darch.learnRateWeights = 0.01,
                 darch.learnRateBiases = 0.01,
                 darch.dropoutInput = 0.,
                 darch.dropoutHidden = 0.,
                 darch.fineTuneFunction = backpropagation, # finetune_SGD,
                 darch.batchSize = 5,
                 darch.numEpochs = 1
)

darch@executeFunction <- runDArch

plot(target, predict(darch))

# run finetune_SGD_bn with batch normalization off

darch@learnRateBiases <- exp(1)
darch@learnRateWeights <- exp(1)

for(i in 1:100) {
  darch <- finetune_SGD_bn(darch, input, target, learn_rate_gamma = exp(-8), with_BN = F)
  ret <- backpropagate_delta_bn(darch, input, target, with_BN = F)
  output <- ret[[4]][[2]]
  delta_weight <- ret[[1]]
  mse_err <- mseError(target, output)
  print(paste0(mse_err[[1]], ": ", mse_err[[2]]))
}

plot(target, output)

darch = darch( x = input,
               y = target,
               layers = c(2, 10, 1),
               # darch = darch,
               darch.layerFunctionDefault = rectified_linear_unit_function,
               darch.layerFunctions = c("2" = linearUnitDerivative),
               darch.bootstrap = F,
               darch.isBin = F,
               darch.isClass = F,
               darch.learnRateWeights = 0.01,
               darch.learnRateBiases = 0.01,
               darch.dropoutInput = 0.,
               darch.dropoutHidden = 0.,
               darch.fineTuneFunction = finetune_SGD_bn, # ,
               with_BN = F,
               darch.batchSize = 5,
               darch.numEpochs = 100
)

darch = darch( x = input,
               y = target,
               layers = c(2, 10, 1),
               # darch = darch,
               darch.layerFunctionDefault = rectified_linear_unit_function,
               darch.layerFunctions = c("2" = linearUnitDerivative),
               darch.bootstrap = F,
               darch.isBin = F,
               darch.isClass = F,
               darch.learnRateWeights = 0.01,
               darch.learnRateBiases = 0.01,
               darch.dropoutInput = 0.,
               darch.dropoutHidden = 0.,
               darch.fineTuneFunction = backpropagation,
               darch.batchSize = 5,
               darch.numEpochs = 100
)


# test run_darch_bn and backpropagate_delta_bn functions

darch = darch( x = input,
               y = target,
               layers = c(2, 10, 1),
               # darch = darch,
               darch.layerFunctionDefault = rectified_linear_unit_function,
               darch.layerFunctions = c("2" = linearUnitDerivative),
               darch.bootstrap = F,
               darch.isBin = F,
               darch.isClass = F,
               darch.learnRateWeights = 0.01,
               darch.learnRateBiases = 0.01,
               darch.dropoutInput = 0.,
               darch.dropoutHidden = 0.,
               darch.fineTuneFunction = finetune_SGD_bn,
               with_BN = T,
               darch.batchSize = 5,
               darch.numEpochs = 1
)

output1 <- predict(darch, newdata = input)
ret <- backpropagate_delta_bn(darch, input, target, with_BN = T)
output2 <- ret[[4]][[2]]
