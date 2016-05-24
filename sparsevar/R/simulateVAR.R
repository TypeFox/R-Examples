#' @title VAR simulation
#'
#' @description This function generates a simulated multivariate VAR time series.
#' 
#' @param N dimension of the time series. 
#' @param p number of lags of the VAR model.
#' @param nobs number of observations to be generated.
#' @param rho base value for the covariance matrix.
#' @param sparsity density (in percentage) of the number of nonzero elements of the VAR matrices.
#' @param method which method to use to generate the VAR matrix. Possible values
#' are \code{"normal"} or \code{"bimodal"}.
#' @param covariance type of covariance matrix to use in the simulation. Possible 
#' values: \code{"toeplitz"}, \code{"block1"}, \code{"block2"} or simply \code{"diagonal"}.
#' 
#' @return A a list of NxN matrices ordered by lag
#' @return data a list with two elements: \code{series} the multivariate time series and
#' \code{noises} the time series of errors
#' @return S the variance/covariance matrix of the process
#'
#' @export
simulateVAR <- function(N = 100, p = 1, nobs = 250, rho = 0.5, sparsity = 0.05, 
                        method = "normal", covariance = "toeplitz") {
  
  # Create the list of the VAR matrices
  A <- list()
  cA <- matrix(0, nrow = N, ncol = N * p)
  for (i in 1:p) {
    A[[i]] <- createSparseMatrix(sparsity = sparsity, N = N, method = method, stationary = TRUE)
    l <- max(Mod(eigen(A[[i]])$values))
    while ((l > 1) | (l == 0)) {
      A[[i]] <- createSparseMatrix(sparsity = sparsity, N = N, method = method, stationary = TRUE)
      l <- max(Mod(eigen(A[[i]])$values))
    }
    A[[i]] <- 1/sqrt(p) * A[[i]]
    cA[1:N, ((i-1) * N) + (1:N)] <- A[[i]]
  }

  # Covariance Matrix: Toeplitz, Block1 or Block2
  if (covariance == "block1"){
    
    l <- floor(N/2)
    I <- diag(1 - rho, nrow = N)
    r <- matrix(0, nrow = N, ncol = N)
    r[1:l, 1:l] <- rho
    T <- I + r
      
  } else if (covariance == "block2") {
  
    l <- floor(N/2)
    I <- diag(1 - rho, nrow = N)
    r <- matrix(0, nrow = N, ncol = N)
    r[1:l, 1:l] <- rho
    r[(l+1):N, (l+1):N] <- rho
    T <- I + r
      
  } else if (covariance == "toeplitz"){
    
    r <- rho^(1:N)
    T <- Matrix::toeplitz(r) 
  
  } else if (covariance == "diagonal"){
    
    T <- diag(x = rho, nrow = N, ncol = N)

  } else {
    
    stop("Unknown covariance matrix type. Possible choices are: toeplitz, block1, block2")
    
  }
  
  # Matrix for MA part
  theta <- matrix(0, N, N)
  ar <- 1:p
  
  # Generate the VAR process 
  data <- MTS::VARMAsim(nobs = nobs, arlags = ar, malags = 0, cnst = 0, phi = cA, theta = theta, skip = 200, sigma = T)
  
  # Output
  out <- list()
  out$data <- data
  out$A <- A
  out$S <- T
  return(out)
  
}
