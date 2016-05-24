#' Calculates the delta functions using backpropagation
#'
#' function that calculates the delta function of a darch object with batch
#' normalization
#'
#' @param darch a darch instance
#' @param trainData training input
#' @param targetData training target
#' @param errorFunc error function to minimize during training. Right now mean squared
#'  erros and cross entropy errors are supported.
#' @param with_BN traing with batch normalization on or off
#'
#' @importFrom darch getLayer
#' @importFrom darch getDropoutMask
#'
#'
#' @references Batch Normalization: Accelerating Deep Network Training by Reducing Internal Covariate Shift
#'  Sergey Ioffe, Christian Szegedy
#' @seealso \url{http://jmlr.org/proceedings/papers/v37/ioffe15.pdf} Pg 4


backpropagate_delta_bn <- function(darch,
                                   trainData,
                                   targetData,
                                   errorFunc = meanSquareErr,
                                   with_BN = TRUE) {

  numLayers <- length(darch@layers)
  layers <- list()
  epsilon <- exp(-12) # a numerical stablizaer used in batch normalization
  numObs <- dim(trainData)[[1]]

  for(i in 1:numLayers) {
    ret <- getLayer(darch, i)[[1]]
    dimV_input <- dim(ret)[[1]] - 1
    dimV_output <- dim(ret)[[2]]

    layers[["weight"]][[i]] <- ret[1:dimV_input, ]

    if(length(getLayer(darch, i)) < 4 | with_BN == FALSE) {
      layers[["gamma"]][[i]] <-
        matrix(rep(1, dimV_output * numObs), numObs, byrow = TRUE)
    } else {
      layers[["gamma"]][[i]] <-
        matrix(rep(getLayer(darch, i)[[4]], numObs), numObs, byrow = TRUE)
    }

    layers[["beta"]][[i]] <- verticalize(ret[(dimV_input + 1),], numObs)

    layers[["x"]][[i]] <- list()
    layers[["mu"]][[i]] <- list()
    layers[["sigma_2"]][[i]] <- list()
    layers[["x_hat"]][[i]] <- list()
    layers[["y"]][[i]] <- list()

    layers[["delta_weight"]][[i]] <- list()
    layers[["delta_x"]][[i]] <- list()
    layers[["delta_y"]][[i]] <- list()
    layers[["delta_beta"]][[i]] <- list()
    layers[["delta_gamma"]][[i]] <- list()
    layers[["output"]][[i]] <- list()
    layers[["derivative"]][[i]] <- list()
  }

  # apply input dropout mask to data
  # TODO same input dropout mask for all data in a batch?
  trainData <- applyDropoutMask(trainData, getDropoutMask(darch, 0))

  # 1. Forwardpropagate
  data <- trainData
  for (i in 1:numLayers){
    weights <- layers[["weight"]][[i]]
    func <- getLayer(darch, i)[[2]]
    # Batch Normalization
    layers[["x"]][[i]] <- data %*% weights

    if(length(getLayer(darch, i)) < 4 | with_BN == FALSE) {
      ret <- batch_normalization(layers[["x"]][[i]],
                                 layers[["gamma"]][[i]],
                                 layers[["beta"]][[i]],
                                 mu = verticalize(rep(0, dim(layers[["gamma"]][[i]])[[2]]), numObs),
                                 sigma_2 = verticalize(rep(1 - epsilon, dim(layers[["gamma"]][[i]])[[2]]), numObs),
                                 epsilon = epsilon
                                 )

    } else {
      ret <- batch_normalization(layers[["x"]][[i]],
                                 layers[["gamma"]][[i]],
                                 layers[["beta"]][[i]],
                                 epsilon = epsilon )
    }
    layers[["mu"]][[i]] <- ret[[1]]
    layers[["sigma_2"]][[i]] <- ret[[2]]
    layers[["x_hat"]][[i]] <- ret[[3]]
    layers[["y"]][[i]] <- ret[[4]]

    ret <- list()

    unit_matrix <- diag(dim(layers[['y']][[i]])[[2]])
    ret <- func(layers[["y"]][[i]],unit_matrix)
    # apply dropout masks to output, unless we're on the last layer
    if (i < numLayers)
    {
      ret[[1]] <- applyDropoutMask(ret[[1]], getDropoutMask(darch, i))
      ret[[2]] <- applyDropoutMask(ret[[2]], getDropoutMask(darch, i))
    }

    layers[["output"]][[i]] <- ret[[1]]
    data <- ret[[1]]
    layers[["derivative"]][[i]] <- ret[[2]]
  }

  # End of forward propagation

  # 2. Calculate the Error on the network output layer
  errorDerivative <- errorFunc(layers[["output"]][[numLayers]], targetData)[[2]]
  layers[["delta_y"]][[numLayers]] <- errorDerivative * layers[["derivative"]][[numLayers]]

  if(length(getLayer(darch, numLayers)) < 4 | with_BN == FALSE) {
    ret <- batch_normalization_differential(layers[["delta_y"]][[numLayers]],
                                            layers[["mu"]][[numLayers]],
                                            layers[["sigma_2"]][[numLayers]],
                                            layers[["x"]][[numLayers]],
                                            layers[["x_hat"]][[numLayers]],
                                            layers[["y"]][[numLayers]],
                                            layers[["gamma"]][[numLayers]],
                                            layers[["beta"]][[numLayers]],
                                            with_BN = FALSE
                                            )

  } else {
    ret <- batch_normalization_differential(layers[["delta_y"]][[numLayers]],
                                            layers[["mu"]][[numLayers]],
                                            layers[["sigma_2"]][[numLayers]],
                                            layers[["x"]][[numLayers]],
                                            layers[["x_hat"]][[numLayers]],
                                            layers[["y"]][[numLayers]],
                                            layers[["gamma"]][[numLayers]],
                                            layers[["beta"]][[numLayers]],
                                            with_BN = TRUE)
  }

  layers[["delta_x"]][[numLayers]] <- ret[[1]]
  layers[["delta_gamma"]][[numLayers]] <- ret[[2]]
  layers[["delta_beta"]][[numLayers]] <- ret[[3]]

  if (numLayers > 1) {
    layers[["delta_weight"]][[numLayers]] <- t(layers[["output"]][[numLayers - 1]]) %*%
      layers[["delta_y"]][[numLayers]]
  } else {
    layers[["delta_weight"]][[numLayers]] <- t(trainData) %*%
      layers[["delta_y"]][[numLayers]]
  }
  # End of calculation

  # 3. Backpropagate the error
  for(i in (numLayers-1):1){
    error <-  layers[["delta_x"]][[i+1]] %*% t(layers[["weight"]][[i + 1]])
    # zero derivatives makes sure that dropout nodes' delta functions are zeros
    layers[["delta_y"]][[i]] <- error * layers[["derivative"]][[i]]

    if(length(getLayer(darch, i)) < 4 | with_BN == FALSE) {
      ret <- batch_normalization_differential(layers[["delta_y"]][[i]],
                                              layers[["mu"]][[i]],
                                              layers[["sigma_2"]][[i]],
                                              layers[["x"]][[i]],
                                              layers[["x_hat"]][[i]],
                                              layers[["y"]][[i]],
                                              layers[["gamma"]][[i]],
                                              layers[["beta"]][[i]],
                                              with_BN = FALSE)

    } else {
      ret <- batch_normalization_differential(layers[["delta_y"]][[i]],
                                              layers[["mu"]][[i]],
                                              layers[["sigma_2"]][[i]],
                                              layers[["x"]][[i]],
                                              layers[["x_hat"]][[i]],
                                              layers[["y"]][[i]],
                                              layers[["gamma"]][[i]],
                                              layers[["beta"]][[i]],
                                              with_BN = TRUE)
    }

    layers[["delta_x"]][[i]] <- ret[[1]]
    layers[["delta_gamma"]][[i]] <- ret[[2]]
    layers[["delta_beta"]][[i]] <- ret[[3]]

    if (i > 1) {
      layers[["delta_weight"]][[i]] <- t(layers[["output"]][[i - 1]]) %*% layers[["delta_y"]][[i]]
    } else {
      layers[["delta_weight"]][[i]] <- t(trainData) %*% layers[["delta_y"]][[i]]
    }

  }

  ret <- list()
  ret[[1]] <- layers[["delta_weight"]]
  ret[[2]] <- layers[["delta_beta"]]
  ret[[3]] <- layers[["delta_gamma"]]
  ret[[4]] <- layers[["output"]]
  ret[[5]] <- layers[["derivative"]]
  ret[[6]] <- layers[["delta_mu"]]
  ret[[7]] <- layers[["delta_sigma_2"]]
  ret[[8]] <- layers[["mu"]]
  ret[[9]] <- layers[["sigma_2"]]
  return(ret)
}
