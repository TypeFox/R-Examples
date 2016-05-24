#' Batch Normalization Function that normalizes the input before applying non-linearity
#'
#' This function normalizes the distribution of inputs to hidden layers in
#' a neural network
#' @param x weighted sum of outputs from the previous layer
#' @param gamma the gamma coefficient
#' @param beta the beta coefficient
#' @param mu the mean of the input neurons. If NULL, it will be caluclated in the function.
#' @param sigma_2 the variance of the input nerurons. If NULL, it will be calcualted in the function.
#' @param epsilon a constant added to the variance for numerical stability
#' @references Batch Normalization: Accelerating Deep Network Training by Reducing Internal Covariate Shift
#' Sergey Ioffe, Christian Szegedy
#' @seealso  \url{http://jmlr.org/proceedings/papers/v37/ioffe15.pdf} Pg 4

batch_normalization <- function(x,
                                gamma,
                                beta,
                                mu = NULL,
                                sigma_2 = NULL,
                                epsilon = exp(-12)) {

  # helper function that repeat a row vector N times
  verticalize <- function(vector, N) {
    return(matrix(rep(vector, N), N, byrow = T))
  }

  numObs <- dim(x)[[1]]
  if(is.null(mu)) {
    mu <-verticalize(colMeans(x), numObs)
  }

  if(is.null(sigma_2)) {
    sigma_2 <- numObs / (numObs - 1) * (verticalize(colMeans(x^2), numObs) - mu^2)
  }



  x_hat <- (x - mu) / sqrt(sigma_2 + epsilon)
  y <- x_hat * gamma + beta

  ret <- list()
  ret[[1]] <- mu
  ret[[2]] <- sigma_2
  ret[[3]] <- x_hat
  ret[[4]] <- y
  return(ret)
}

#' Function that calcualtes the differentials in the batch normalization mode
#'
#' Calculates the differentials in batch normalization
#'
#' @param delta_y derivative wrt y
#' @param mu mean of the input
#' @param sigma_2 variance of the input
#' @param x input
#' @param x_hat normalized input
#' @param y transformed input after batch normalization
#' @param gamma gamma coefficient
#' @param beta beta coefficient
#' @param epsilon the contant added to the variance for numeric stability
#' @param with_BN logical value, set to TRUE to turn on batch normalization
#'
#' @references Batch Normalization: Accelerating Deep Network Training by Reducing Internal Covariate Shift
#'  Sergey Ioffe, Christian Szegedy
#' @seealso \url{http://jmlr.org/proceedings/papers/v37/ioffe15.pdf} Pg 4

batch_normalization_differential <- function(delta_y,
                                             mu,
                                             sigma_2,
                                             x,
                                             x_hat,
                                             y,
                                             gamma,
                                             beta,
                                             epsilon = exp(-12),
                                             with_BN = T) {
  # helper function that repeat a row vector N times
  verticalize <- function(vector, N) {
    return(matrix(rep(vector, N), N, byrow = T))
  }
  numObs <- dim(x)[[1]]

  delta_x_hat <- delta_y * gamma

  if(with_BN) {
    delta_sigma_2 <-  verticalize(colSums(delta_x_hat * (x - mu) * (-0.5) * (sigma_2 + epsilon)^(-1.5)), numObs)

    tmp1 <- verticalize(colSums(delta_x_hat * (-1) / sqrt(sigma_2 + epsilon)), numObs)
    tmp2 <- delta_sigma_2 * verticalize(colMeans(-2 * (x- mu)), numObs)

    delta_mu <- tmp1 + tmp2

    delta_gamma <- verticalize(colSums(delta_y * x_hat), numObs)
  } else {
    delta_sigma_2 <-  verticalize(rep(0, dim(delta_y)[[2]]), numObs)
    delta_mu <- verticalize(rep(0, dim(delta_y)[[2]]), numObs)
    delta_gamma <- verticalize(rep(0, dim(delta_y)[[2]]), numObs)
  }

  tmp1 <- delta_x_hat / sqrt(sigma_2 + epsilon)
  tmp2 <- delta_sigma_2 * 2 * (x - mu) / numObs
  tmp3 <- delta_mu / numObs
  delta_x <- tmp1 + tmp2 + tmp3

  delta_beta <- verticalize(colSums(delta_y), numObs)

  ret <- list()
  ret[[1]] <- delta_x
  ret[[2]] <- delta_gamma
  ret[[3]] <- delta_beta
  ret[[4]] <- delta_x_hat
  ret[[5]] <- delta_sigma_2
  ret[[6]] <- delta_mu
  return(ret)
}
