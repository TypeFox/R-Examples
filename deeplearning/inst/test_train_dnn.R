###########################################################################################################################
# Test case 1: test the basic functionality of run_dnn

input <- matrix(runif(100), 50, 2)
input_valid <- matrix(runif(10), 5, 2)
target <- rowSums(input + input^2)
target_valid <- rowSums(input_valid + input_valid^2)

darch <- new_dnn(c(2,5,1))
darch <- train_dnn(darch,
          input,
          target,
          # input_valid,
          # target_valid,

          learn_rate_weight = exp(-5),
          learn_rate_bias = exp(-5),
          learn_rate_gamma = exp(-5),
          batch_size = 10,
          batch_normalization = F,
          dropout_input = 0,
          dropout_hidden = 0,
          momentunm_initial = 0.6,
          momentum_final = 0.9,
          momentum_switch = 100,
          num_epochs = 100,

          # target types
          error_function = meanSquareErr,
          report_classification_error = F
)

# test the dropout


###########################################################################################################################
# Test case 2: Test mixed training of BN and No BN
# 2.1
# first train with BN on
# then train with BN off
rm(darch)
darch <- new_dnn(c(2, 5, 10, 1))
darch <- train_dnn(darch,
                   input,
                   target,
                   input_valid,
                   target_valid,
                   # training parameters
                   learn_rate_weight = exp(-5),
                   learn_rate_bias = exp(-5),
                   learn_rate_gamma = exp(-5),
                   batch_size = 10,
                   batch_normalization = T,
                   dropout_input = 0.4,
                   dropout_hidden = 0.8,
                   momentunm_initial = 0.6,
                   momentum_final = 0.9,
                   momentum_switch = 100,
                   num_epochs = 50,
                   # target types
                   error_function = meanSquareErr,
                   report_classification_error = F
)

darch <- train_dnn(darch,
                   input,
                   target,
                   input_valid,
                   target_valid,
                   # training parameters
                   learn_rate_weight = exp(-5),
                   learn_rate_bias = exp(-5),
                   learn_rate_gamma = exp(-5),
                   batch_size = 10,
                   batch_normalization = F,
                   dropout_input = 0.4,
                   dropout_hidden = 0.8,
                   momentunm_initial = 0.6,
                   momentum_final = 0.9,
                   momentum_switch = 100,
                   num_epochs = 50,
                   # target types
                   error_function = meanSquareErr,
                   report_classification_error = F
)

plot(darch@stats$dataErrors$raw)
plot(darch@stats$validErrors$raw)

# 2.2
# firt train with BN off
# then trian with BN on

rm(darch)
darch <- new_dnn(c(2, 5, 10, 1))
darch <- train_dnn(darch,
                   input,
                   target,
                   input_valid,
                   target_valid,
                   # training parameters
                   learn_rate_weight = exp(-5),
                   learn_rate_bias = exp(-5),
                   learn_rate_gamma = exp(-5),
                   batch_size = 10,
                   batch_normalization = F,
                   dropout_input = 0.4,
                   dropout_hidden = 0.8,
                   momentunm_initial = 0.6,
                   momentum_final = 0.9,
                   momentum_switch = 100,
                   num_epochs = 50,
                   # target types
                   error_function = meanSquareErr,
                   report_classification_error = F
)


darch <- train_dnn(darch,
                   input,
                   target,
                   input_valid,
                   target_valid,
                   # training parameters
                   learn_rate_weight = exp(-5),
                   learn_rate_bias = exp(-5),
                   learn_rate_gamma = exp(-5),
                   batch_size = 10,
                   batch_normalization = T,
                   dropout_input = 0.4,
                   dropout_hidden = 0.8,
                   momentunm_initial = 0.6,
                   momentum_final = 0.9,
                   momentum_switch = 100,
                   num_epochs = 50,
                   # target types
                   error_function = meanSquareErr,
                   report_classification_error = F
)

plot(darch@stats$dataErrors$raw)
plot(darch@stats$validErrors$raw)


###########################################################################################################################
# Test 3: Evaluate the Batch Normalization
# Compare BN training with no BN training

input <- matrix(runif(500), 250, 2)
input_valid <- matrix(runif(50), 25, 2)
target <- rowSums(cos(input) + sin(input)^2 + tan(input)^3)
target_valid <- as.matrix(rowSums(cos(input_valid) + sin(input_valid)^2 + tan(input_valid)^3))


rm(darch)
darch <- new_dnn(c(2, 10, 10, 1))
darch <- train_dnn(darch,
                   input,
                   target,
                   input_valid,
                   target_valid,
                   # training parameters
                   learn_rate_weight = exp(-8),
                   learn_rate_bias = exp(-8),
                   learn_rate_gamma = exp(-8),
                   batch_size = 10,
                   batch_normalization = T,
                   dropout_input = 0.,
                   dropout_hidden = 0.,
                   momentunm_initial = 0.6,
                   momentum_final = 0.9,
                   momentum_switch = 100,
                   num_epochs = 250,
                   # target types
                   error_function = meanSquareErr,
                   report_classification_error = F
)

rsq(darch, input = input, target = target)
# 50 Iter: .959
# 250 Iter: .961
lines(x = c(min(target), max(target)), y =  c(min(target), max(target)), col = "red")

rsq(darch, input = input_valid, target = target_valid)
# 50 iterations: .965
# 250 iterations: .964
lines(x = c(min(target), max(target)), y =  c(min(target), max(target)), col = "red")

plot(darch@stats$dataErrors$raw)

plot(darch@stats$validErrors$raw)


rm(darch_2)
darch_2 <- new_dnn(c(2, 10, 10, 1))
darch_2 <- train_dnn(darch_2,
                   input,
                   target,
                   input_valid,
                   target_valid,
                   # training parameters
                   learn_rate_weight = exp(-8),
                   learn_rate_bias = exp(-8),
                   learn_rate_gamma = exp(-8),
                   batch_size = 10,
                   batch_normalization = F,
                   dropout_input = 0.,
                   dropout_hidden = 0.,
                   momentunm_initial = 0.6,
                   momentum_final = 0.9,
                   momentum_switch = 100,
                   num_epochs = 250,
                   # target types
                   error_function = meanSquareErr,
                   report_classification_error = F
)

rsq(darch_2, input = input, target = target)
# 50 iterations:  .687
# 250 iterations: .780
lines(x = c(min(target), max(target)), y =  c(min(target), max(target)), col = "red")

rsq(darch_2, input = input_valid, target = target_valid)
# 50 iterations: .728
# 250 iterations: .881
lines(x = c(min(target), max(target)), y =  c(min(target), max(target)), col = "red")

plot(darch_2@stats$dataErrors$raw)

plot(darch_2@stats$validErrors$raw)



###########################################################################################################################
# Test case 4:
# Cross comparison between train_dnn and darch
# differences:
# 1) batch generation
# 2) batch normalization
# 3) Cross Entropy Error
# 4) Bug in runDArch with dropout

# 4.1 benchmark - 1 batch. use no batch normalization

input <- matrix(runif(500), 250, 2)
input_valid <- matrix(runif(50), 25, 2)
target <- rowSums(cos(input) + sin(input)^2)
target_valid <- as.matrix(rowSums(cos(input_valid) + sin(input_valid)^2))

# use train_dnn function from deeplearning library

rm(darch)
darch <- new_dnn(c(2, 20, 20, 1))
darch <- train_dnn(darch,
                   input,
                   target,
                   input_valid,
                   target_valid,
                   # training parameters
                   learn_rate_weight = exp(-8),
                   learn_rate_bias = exp(-8),
                   learn_rate_gamma = exp(-8),
                   batch_size = 250,
                   batch_normalization = F,
                   dropout_input = 0.,
                   dropout_hidden = 0.,
                   momentunm_initial = 0.6,
                   momentum_final = 0.9,
                   momentum_switch = 100,
                   num_epochs = 500,
                   # target types
                   error_function = meanSquareErr,
                   report_classification_error = F
)


rsq(darch, input = input_valid, target = target_valid) # .983
lines(x = c(2,3), y = c(2, 3), col = "red")

plot(darch@stats$dataErrors$raw)

plot(darch@stats$validErrors$raw)

# use darch function from the darch library

rm(darch)

darch <- darch(    input,
                   target,
                   layers = c(2, 20, 20, 1),
                   xvalid = input_valid,
                   yValid = target_valid,
                   # training parameters
                   darch.learnRateBiases = exp(-8),
                   darch.learnRateWeights = exp(-8),
                   darch.layerFunctionDefault = rectified_linear_unit_function,
                   darch.layerFunctions = list("3" = linearUnitDerivative),
                   darch.batchSize = 250,
                   darch.dropoutInput = 0.,
                   darch.dropoutHidden = 0.,
                   darch.momentumSwitch = 100,
                   darch.initialMomentum = 0.6,
                   darch.finalMomentum = 0.9,
                   darch.numEpochs = 500,
                   darch.isBin = F,
                   darch.isClass = F
)

rsq(darch, input = input_valid, target = target_valid) # .986
lines(x = c(2,3), y = c(2, 3), col = "red")

plot(darch@stats$dataErrors$raw)

plot(darch@stats$validErrors$raw)


########################################################################
# 4.2 test batch initialization - 50 batches. use no batch normalization


input <- matrix(runif(500), 250, 2)
input_valid <- matrix(runif(50), 25, 2)
target <- rowSums(cos(input) + sin(input)^2 + tan(input)^3)
target_valid <- as.matrix(rowSums(cos(input_valid) + sin(input_valid)^2 + tan(input_valid)^3))

# use train_dnn function from deeplearning library

rm(darch)
darch <- new_dnn(c(2, 20, 20, 1))
darch <- train_dnn(darch,
                   input,
                   target,
                   input_valid,
                   target_valid,
                   # training parameters
                   learn_rate_weight = exp(-8),
                   learn_rate_bias = exp(-8),
                   learn_rate_gamma = exp(-8),
                   batch_size = 10,
                   batch_normalization = F,
                   dropout_input = 0.,
                   dropout_hidden = 0.,
                   momentunm_initial = 0.6,
                   momentum_final = 0.9,
                   momentum_switch = 100,
                   num_epochs = 500,
                   # target types
                   error_function = meanSquareErr,
                   report_classification_error = F
)

rsq(darch, input = input, target = target)
# 100 iterations: .760
# 500 iterations: .987
lines(x = c(min(target), max(target)), y =  c(min(target), max(target)), col = "red")

rsq(darch, input = input_valid, target = target_valid)
# 100 iterations: .770
# 500 iterations; .979
lines(x = c(min(target), max(target)), y =  c(min(target), max(target)), col = "red")

plot(darch@stats$dataErrors$raw)

plot(darch@stats$validErrors$raw)

# use darch function from the darch library

darch_2 <- darch(    input,
                     target,
                     layers = c(2, 20, 20, 1),
                     xvalid = input_valid,
                     yValid = target_valid,
                     # training parameters
                     darch.learnRateBiases = exp(-8),
                     darch.learnRateWeights = exp(-8),
                     darch.layerFunctionDefault = rectified_linear_unit_function,
                     darch.layerFunctions = list("3" = linearUnitDerivative),
                     darch.batchSize = 10,
                     darch.dropoutInput = 0.,
                     darch.dropoutHidden = 0.,
                     darch.momentumSwitch = 100,
                     darch.initialMomentum = 0.6,
                     darch.finalMomentum = 0.9,
                     darch.numEpochs = 500,
                     darch.isBin = F,
                     darch.isClass = F
)

rsq(darch_2)
# 100 iterations: .767
# 500 iterations: .980
lines(x = c(min(target), max(target)), y =  c(min(target), max(target)), col = "red")

rsq(darch_2, input = input_valid, target = target_valid)
# 100 iterations: .733
# 500 iterations: .974
lines(x = c(min(target), max(target)), y =  c(min(target), max(target)), col = "red")

plot(darch_2@stats$dataErrors$raw)

plot(darch_2@stats$validErrors$raw)



####################################################################
# 4.3 test batchnormalization - 50 batches. use  batch normalization

input <- matrix(runif(500), 250, 2)
input_valid <- matrix(runif(50), 25, 2)
target <- rowSums(cos(input) + sin(input)^2 + tan(input)^3)
target_valid <- as.matrix(rowSums(cos(input_valid) + sin(input_valid)^2 + tan(input_valid)^3))

# use train_dnn function from deeplearning library

rm(darch)
darch <- new_dnn(c(2, 20, 20, 1))
darch <- train_dnn(darch,
                   input,
                   target,
                   input_valid,
                   target_valid,
                   # training parameters
                   learn_rate_weight = exp(-8),
                   learn_rate_bias = exp(-8),
                   learn_rate_gamma = exp(-8),
                   batch_size = 10,
                   batch_normalization = T,
                   dropout_input = 0.,
                   dropout_hidden = 0,
                   momentunm_initial = 0.6,
                   momentum_final = 0.9,
                   momentum_switch = 100,
                   num_epochs = 100,
                   # target types
                   error_function = meanSquareErr,
                   report_classification_error = F
)

rsq(darch, input = input, target = target)
# 100 Iterations: .968
# 500 Iterations: .971
lines(x = c(min(target), max(target)), y =  c(min(target), max(target)), col = "red")

rsq(darch, input = input_valid, target = target_valid)
# 100 Iterations: .946
# 500 Iterations: .930
lines(x = c(min(target), max(target)), y =  c(min(target), max(target)), col = "red")

plot(darch@stats$dataErrors$raw)

plot(darch@stats$validErrors$raw)




# use darch function from the darch library
darch_2 <- darch(    input,
                     target,
                     layers = c(2, 20, 20, 1),
                     xvalid = input_valid,
                     yValid = target_valid,
                     # training parameters
                     darch.learnRateBiases = exp(-8),
                     darch.learnRateWeights = exp(-8),
                     darch.layerFunctionDefault = rectified_linear_unit_function,
                     darch.layerFunctions = list("3" = linearUnitDerivative),
                     darch.batchSize = 10,
                     darch.dropoutInput = 0.,
                     darch.dropoutHidden = 0,
                     darch.momentumSwitch = 100,
                     darch.initialMomentum = 0.6,
                     darch.finalMomentum = 0.9,
                     darch.numEpochs = 500,
                     darch.isBin = F,
                     darch.isClass = F
)

rsq(darch_2)
# 100 Iterations: .727
# 500 Iterations: .974
lines(x = c(min(target), max(target)), y =  c(min(target), max(target)), col = "red")

rsq(darch_2, input = input_valid, target = target_valid)
# 100 Iterations: .742
# 500 Iterations: .974
lines(x = c(min(target), max(target)), y =  c(min(target), max(target)), col = "red")

plot(darch_2@stats$dataErrors$raw)

plot(darch_2@stats$validErrors$raw)



#################################################################################
# 4.4 test batchnormalization - 50 batches. use  batch normalization. use dropout

input <- matrix(runif(500), 250, 2)
input_valid <- matrix(runif(50), 25, 2)
target <- rowSums(cos(input) + sin(input)^2 + tan(input)^3)
target_valid <- as.matrix(rowSums(cos(input_valid) + sin(input_valid)^2 + tan(input_valid)^3))

# use train_dnn function from deeplearning library

rm(darch)
darch <- new_dnn(c(2, 40, 40, 1), hidden_layer_default = sigmoidUnitDerivative)
darch <- train_dnn(darch,
                   input,
                   target,
                   input_valid,
                   target_valid,
                   # training parameters

                   learn_rate_weight = exp(-8) * 100,
                   learn_rate_bias = exp(-8) * 100,
                   learn_rate_gamma = exp(-8) * 100,
                   batch_size = 10,
                   batch_normalization = T,
                   dropout_input = 0.2,
                   dropout_hidden = 0,
                   momentunm_initial = 0.6,
                   momentum_final = 0.9,
                   momentum_switch = 100,
                   num_epochs = 500,
                   # target types
                   error_function = meanSquareErr,
                   report_classification_error = F
)

rsq(darch, input = input, target = target)
# learn rate: exp(-8) * 100
# dropout input/hidden: .2/.3
# 100 Iterations: .742
# 500 Iterations: .937
lines(x = c(min(target), max(target)), y =  c(min(target), max(target)), col = "red")

rsq(darch, input = input_valid, target = target_valid)
# learn rate: exp(-8) * 100
# dropout input/hidden: .2/.3
# 100 Iterations: .782
# 500 Iterations: .914
lines(x = c(min(target), max(target)), y =  c(min(target), max(target)), col = "red")

plot(darch@stats$dataErrors$raw)

plot(darch@stats$validErrors$raw)

# use darch function from the darch library



darch_2 <- darch(    input,
                     target,
                     layers = c(2, 40, 40, 1),
                     xvalid = input_valid,
                     yValid = target_valid,
                     # training parameters
                     darch.learnRateBiases = exp(-8) * 100,
                     darch.learnRateWeights = exp(-8) * 100,
                     darch.layerFunctionDefault = sigmoidUnitDerivative,
                     darch.layerFunctions = list("3" = linearUnitDerivative),
                     darch.batchSize = 10,
                     darch.dropoutInput = 0.2,
                     darch.dropoutHidden = 0.3,
                     darch.momentumSwitch = 100,
                     darch.initialMomentum = 0.6,
                     darch.finalMomentum = 0.9,
                     darch.numEpochs = 100,
                     darch.isBin = F,
                     darch.isClass = F
)

# drop out fails with ReLU!!!!!!!!!!!!!!!!!!!!!!

rsq(darch_2, input = input, target = target)
# learn rate: exp(-8) * 100
# dropout input/hidden: .2/.3
# 100 Iterations: .638
# 500 Iterations: .657
lines(x = c(min(target), max(target)), y =  c(min(target), max(target)), col = "red")

rsq(darch_2, input = input_valid, target = target_valid)
# learn rate: exp(-8) * 100
# dropout input/hidden: .2/.3
# 100 Iterations: .600
# 500 Iterations: .347
lines(x = c(min(target), max(target)), y =  c(min(target), max(target)), col = "red")

plot(darch_2@stats$dataErrors$raw)

plot(darch_2@stats$validErrors$raw)




#################################################################################
# 4.5 test error functions

input <- matrix(runif(500), 250, 2)
input_valid <- matrix(runif(50), 25, 2)
target <- rowSums(cos(input) + sin(input)^2 + tan(input)^3)
target_valid <- as.matrix(rowSums(cos(input_valid) + sin(input_valid)^2 + tan(input_valid)^3))
med <- median(target)
target <- 1 * (target < (med * runif(1) * 2 ))
target_valid <- 1 * (target_valid < (med * runif(1) * 2))


# use train_dnn function from deeplearning library
rm(darch)
darch <- new_dnn(c(2, 20, 20, 1), output_layer_default = sigmoidUnitDerivative)
darch <- train_dnn(darch,
                   input,
                   target,
                   input_valid,
                   target_valid,
                   # training parameters

                   learn_rate_weight = exp(-8) ,
                   learn_rate_bias = exp(-8),
                   learn_rate_gamma = exp(-8),
                   batch_size = 10,
                   batch_normalization = T,
                   dropout_input = 0.,
                   dropout_hidden = 0.,
                   momentunm_initial = 0.6,
                   momentum_final = 0.9,
                   momentum_switch = 100,
                   num_epochs = 50,
                   # target types
                   error_function = crossEntropyErr,
                   report_classification_error = T
)

AR(darch, input = input, target = target)

# 100 Iterations:  .916
# 500 Iterations:

AR(darch, input = input_valid, target = target_valid)

# 100 Iterations:  1
# 500 Iterations:

plot(darch@stats$dataErrors$raw)

plot(darch@stats$validErrors$raw)

# use darch function from the darch library



darch_2 <- darch(    input,
                     target,
                     layers = c(2, 20, 20, 1),
                     xvalid = input_valid,
                     yValid = target_valid,
                     # training parameters
                     darch.learnRateBiases = exp(-8) * 1,
                     darch.learnRateWeights = exp(-8) * 1,
                     darch.layerFunctionDefault = rectified_linear_unit_function,
                     darch.layerFunctions = list("3" = sigmoidUnitDerivative),
                     darch.batchSize = 10,
                     darch.dropoutInput = 0.,
                     darch.dropoutHidden = 0.,
                     darch.momentumSwitch = 100,
                     darch.initialMomentum = 0.6,
                     darch.finalMomentum = 0.9,
                     darch.numEpochs = 100,
                     darch.isBin = T,
                     darch.isClass = F
)

# drop out fails with ReLU!!!!!!!!!!!!!!!!!!!!!!

AR(darch_2, input = input, target = target)

# 100 Iterations: .94
# 500 Iterations:

AR(darch_2, input = input_valid, target = target_valid)

# 100 Iterations:  1
# 500 Iterations:

plot(darch_2@stats$dataErrors$raw)

plot(darch_2@stats$validErrors$raw)





