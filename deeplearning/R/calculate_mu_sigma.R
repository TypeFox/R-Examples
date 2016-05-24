#' Calculates the mu and sigmas of a darch instance
#'
#' This function calculates the mu and sigmas of hidden layers in a darch instance
#' @param darch a darch instance
#' @param input input data
#'
#' @importFrom darch getLayer
#'
#'




calcualte_population_mu_sigma <- function (darch, input) {
  numLayers <- length(darch@layers)
  layers <- list()
  epsilon <- exp(-12) # a numerical stablizaer used in batch normalization
  numObs <- dim(input)[[1]]

  for(i in 1:numLayers) {
    ret <- getLayer(darch, i)[[1]]
    dimV_input <- dim(ret)[[1]] - 1
    dimV_output <- dim(ret)[[2]]

    layers[["weight"]][[i]] <- ret[1:dimV_input, ]

    layers[["gamma"]][[i]] <-
        matrix(rep(getLayer(darch, i)[[4]], numObs), numObs, byrow = T)

    layers[["beta"]][[i]] <- verticalize(ret[(dimV_input + 1),], numObs)

    layers[["x"]][[i]] <- list()
    layers[["mu"]][[i]] <- list()
    layers[["sigma_2"]][[i]] <- list()
    layers[["x_hat"]][[i]] <- list()
    layers[["y"]][[i]] <- list()
  }

  # Forwardpropagate
  data <- input
  for (i in 1:numLayers){
    weights <- layers[["weight"]][[i]]
    func <- getLayer(darch, i)[[2]]
    # Batch Normalization
    layers[["x"]][[i]] <- data %*% weights

    ret <- batch_normalization(layers[["x"]][[i]],
                                 layers[["gamma"]][[i]],
                                 layers[["beta"]][[i]],
                                 epsilon = epsilon )

    layers[["mu"]][[i]] <- ret[[1]]
    layers[["sigma_2"]][[i]] <- ret[[2]]
    layers[["x_hat"]][[i]] <- ret[[3]]
    layers[["y"]][[i]] <- ret[[4]]

    ret <- list()

    unit_matrix <- diag(dim(layers[['y']][[i]])[[2]])
    ret <- func(layers[["y"]][[i]],unit_matrix)

    layers[["output"]][[i]] <- ret[[1]]
    data <- ret[[1]]
    layers[["derivative"]][[i]] <- ret[[2]]
  }

  # End of forward propagation

  for (i in 1:numLayers) {
    darch@layers[[i]][[5]] <- layers[["mu"]][[i]][1, ]
    darch@layers[[i]][[6]] <- layers[["sigma_2"]][[i]][1, ]
  }

  return (darch)
}

#' Resets the mu and sigmas of a darch instance to 0 and 1
#'
#' This function resets the mu and sigmas of hidden layers in a darch instance
#'  to 0 and 1
#' @param darch a darch instance
#'
#' @importFrom darch getLayer
#'



reset_population_mu_sigma <- function (darch) {
  numLayers <- length(darch@layers)
  epsilon <- exp(-12) # a numerical stablizaer used in batch normalization

  for(i in 1:numLayers) {
    ret <- getLayer(darch, i)[[1]]
    dimV_output <- dim(ret)[[2]]
    darch@layers[[i]][[5]] <- rep(0, dimV_output)
    darch@layers[[i]][[6]] <- rep(1 - epsilon, dimV_output)
  }

  return (darch)
}

