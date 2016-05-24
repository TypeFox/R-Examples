#' Execution function that runs in the batch normalization mode
#'
#' This function calcualtes the output of a deep neural network with input data
#'
#' @param darch a darch instance
#' @param data input data


run_dnn <- function(darch, data){
  darch@executeOutput <- list()
  layers <- darch@layers
  # If there's only one row of input data, convert vector to matrix
  # TODO make sure that data is matrix before passing it to this function
  if(is.null(dim(data))){
    data <- t(as.matrix(data))
  }

  numRows <- dim(data)[1]

  output <- list()
  derivative <- list()

  for(i in 1:length(layers)){
    ret <- layers[[i]][[1]]
    dimV_input <- dim(ret)[[1]] - 1
    dimV_output <- dim(ret)[[2]]

    weight <- ret[1:(dimV_input), ]
    beta <- verticalize(ret[(dimV_input + 1), ], numRows)

    gamma <- darch@layers[[i]][[4]]
    gamma <- verticalize(gamma, numRows)

    x <- data %*% weight

    mu <- verticalize(layers[[i]][[5]], numRows)

    sigma_2 <- verticalize(layers[[i]][[6]], numRows)

    ret <- batch_normalization(x, gamma, beta, mu, sigma_2)

    y <- ret[[4]]

    unit_matrix <- diag(dim(y)[[2]])
    ret <- layers[[i]][[2]](y, unit_matrix)
    data <- ret[[1]]
    output[[i]] <- ret[[1]]
    derivative[[i]] <- ret[[2]]
  }

  darch@executeOutput <- output
  return(darch)
}
