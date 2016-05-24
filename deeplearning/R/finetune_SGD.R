#' Updates a deep neural network's parameters using stochastic gradient descent
#'  method and batch normalization
#'
#' This function finetunes a DArch network using SGD approach
#'
#' @param darch a darch instance
#' @param trainData training input
#' @param targetData training target
#' @param learn_rate_weight leanring rate for the weight matrices
#' @param learn_rate_bias learning rate for the biases
#' @param learn_rate_gamma learning rate for the gammas
#' @param errorFunc the error function to minimize during training
#' @param with_BN logical value, T to train the neural net with batch normalization
#'
#' @importFrom darch getLayer getDropoutMask getMomentum
#'
#' @return a darch instance with parameters updated with stochastic gradient descent
#'

finetune_SGD_bn <- function(darch,
                            trainData,
                            targetData,
                            learn_rate_weight = exp(-10),
                            learn_rate_bias = exp(-10),
                            learn_rate_gamma = exp(-10),
                            errorFunc = meanSquareErr,
                            with_BN = T) {
    # stats <- getStats(darch)

    ret <- backpropagate_delta_bn(darch, trainData, targetData, errorFunc, with_BN)
    delta_weight <- ret[[1]]
    delta_beta <- ret[[2]]
    delta_gamma <- ret[[3]]

    learnRateBiases <- learn_rate_bias
    learnRateWeights <- learn_rate_weight
    learnRateGamma <- learn_rate_gamma

    numLayers <- length(delta_weight)

    for(i in numLayers:1) {
      weights <- getLayer(darch, i)[[1]]
      biases <- weights[nrow(weights),,drop=F]
      weights <- weights[1:(nrow(weights)-1),,drop=F]
      gamma <- getLayer(darch, i)[[4]]
      weightsChange_prev <- getLayer(darch, i)[[3]]

    # Calculate the change in weights
    # apply dropout mask to momentum
      weightsInc <- (learnRateWeights * delta_weight[[i]])
      weightsChange <- weightsInc + (getMomentum(darch) *
                             weightsChange_prev * getDropoutMask(darch, i-1)
      )
      weights <- weights - weightsChange

      # Calculate the change in beta (biases)
      biasesInc <- learnRateBiases * delta_beta[[i]][1,]
      biases <- biases - biasesInc

      # Calculate the change in gamma
      gammaInc <- learnRateGamma * delta_gamma[[i]][1,]
      gamma <- gamma - gammaInc

      darch@layers[[i]][[1]] <- rbind(weights,biases)
      darch@layers[[i]][[3]] <- weightsInc
      darch@layers[[i]][[4]] <- gamma
  }

  # setStats(darch) <- stats
  return(darch)

}
