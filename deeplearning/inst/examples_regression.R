# Examples of Regression

input <- matrix(runif(1000), 500, 2)
input_valid <- matrix(runif(100), 50, 2)
target <- rowSums(input + input^2)
target_valid <- rowSums(input_valid + input_valid^2)


# create a new deep neural network for classificaiton
dnn_regression <- new_dnn(
                          c(2, 50, 50, 20, 1),  # The layer structure of the deep neural network.
                                                # The first element is the number of input variables.
                                                # The last element is the number of output variables.
                          hidden_layer_default = rectified_linear_unit_function, # for hidden layers, use rectified_linear_unit_function
                          output_layer_default = linearUnitDerivative # for regression, use linearUnitDerivative function
                          )

dnn_regression <- train_dnn(
                     dnn_regression,

                     # training data
                     input, # input variable for training
                     target, # target variable for training
                     input_valid, # input variable for validation
                     target_valid, # target variable for validation

                     # training parameters
                     learn_rate_weight = exp(-8) * 1, # learning rate for weights, higher if use dropout
                     learn_rate_bias = exp(-8) * 1, # learning rate for biases, hihger if use dropout
                     learn_rate_gamma = exp(-8) * 1, # learning rate for the gamma factor used
                     batch_size = 10, # number of observations in a batch during training. Higher for faster training. Lower for faster convergence
                     batch_normalization = T, # logical value, T to use batch normalization
                     dropout_input = 0.2, # dropout ratio in input.
                     dropout_hidden = 0.5, # dropout ratio in hidden layers
                     momentum_initial = 0.6, # initial momentum in Stochastic Gradient Descent training
                     momentum_final = 0.9, # final momentum in Stochastic Gradient Descent training
                     momentum_switch = 100, # after which the momentum is switched from initial to final momentum
                     num_epochs = 100, # number of iterations in training

                     # Error function
                     error_function = meanSquareErr, # error function to minimize during training. For regression, use meanSquareErr
                     report_classification_error = F # whether to print classification error during training
)

# the prediciton by dnn_regression
pred <- predict(dnn_regression)

# calculate the r-squared of the prediciton
rsq(dnn_regression)

# calcualte the r-squared of the prediciton in validation
rsq(dnn_regression, input = input_valid, target = target_valid)

# print the layer weights
# this function can print heatmap, histogram, or a surface
print_weight(dnn_regression, 1, type = "heatmap")

print_weight(dnn_regression, 2, type = "surface")

print_weight(dnn_regression, 3, type = "histogram")




