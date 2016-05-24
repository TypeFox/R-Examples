#' Computes the maximum likelihood estimator for the sample covariance matrix
#' under the assumption of multivariate normality.
#'
#' For a sample matrix, \code{x}, we compute the sample covariance matrix of the
#' data as the maximum likelihood estimator (MLE) of the population covariance
#' matrix.
#'
#' If the \code{diag} option is set to \code{TRUE}, then we assume the population
#' covariance matrix is diagonal, and the MLE is computed under this assumption.
#' In this case, we return a vector of length \code{p} instead.
#'
#' @export
#' @param x data matrix with \code{n} observations and \code{p} feature vectors
#' @param diag logical value. If TRUE, assumes the population covariance matrix
#' is diagonal. By default, we assume that \code{diag} is \code{FALSE}.
#' @return sample covariance matrix of size \eqn{p \times p}. If \code{diag} is
#' \code{TRUE}, then a vector of length \code{p} is returned instead.
cov_mle <- function(x, diag = FALSE) {
  x <- as.matrix(x)
  n <- nrow(x)
  if (diag) {
    (n - 1) / n * apply(x, 2, var)
  } else {
    (n - 1) / n * cov(x)
  }
}

#' Computes the pooled maximum likelihood estimator (MLE) for the common
#' covariance matrix
#'
#' For the matrix \code{x}, we compute the MLE for the population covariance
#' matrix under the assumption that the data are sampled from \eqn{K}
#' multivariate normal populations having equal covariance matrices.
#'
#' @export
#' @param x data matrix with \code{n} observations and \code{p} feature vectors
#' @param y class labels for observations (rows) in \code{x}
#' @return pooled sample covariance matrix of size \eqn{p \times p}
#' @examples
#' cov_pool(iris[, -5], iris$Species)
cov_pool <- function(x, y) {
  x <- as.matrix(x)
  y <- as.factor(y)
  n <- length(y)
  
  scatter_matrices <- tapply(seq_len(n), y, function(i) {
    (length(i) - 1) * cov(as.matrix(x[i, ]))
  })
  as.matrix(Reduce("+", scatter_matrices) / n)
}

#' Computes the covariance-matrix maximum likelihood estimators for each class
#' and returns a list.
#'
#' For a sample matrix, \code{x}, we compute the MLE for the covariance matrix
#' for each class given in the vector, \code{y}.
#' @export
#' @param x data matrix with \code{n} observations and \code{p} feature vectors
#' @param y class labels for observations (rows) in \code{x}
#' @return list of the sample covariance matrices of size \eqn{p \times p} for
#' each class given in \code{y}.
cov_list <- function(x, y) {
  x <- as.matrix(x)
  y <- as.factor(y)
  tapply(seq_along(y), y, function(i) {
    cov_mle(x[i, ])
  })
}

#' Computes the eigenvalue decomposition of the maximum likelihood estimators
#' (MLE) of the covariance matrices for the given data matrix
#'
#' For the classes given in the vector \code{y}, we compute the eigenvalue
#' (spectral) decomposition of the class sample covariance matrices (MLEs) using
#' the data matrix \code{x}.
#'
#' If the \code{fast} argument is selected, we utilize the so-called Fast
#' Singular Value Decomposition (SVD) to quickly compute the eigenvalue
#' decomposition. To compute the Fast SVD, we use the \code{\link{fast.svd}}
#' function, which employs a well-known trick for tall data (large \code{n},
#' small \code{p}) and wide data (large \code{p}, small \code{n}) to compute the
#' SVD corresponding to the nonzero singular values. For more information about
#' the Fast SVD, see \code{\link[corpcor]{fast.svd}}.
#' 
#' @importFrom corpcor fast.svd
#' @export
#' @param x data matrix with \code{n} observations and \code{p} feature vectors
#' @param y class labels for observations (rows) in \code{x}
#' @param pool logical. Should the sample covariance matrices be pooled?
#' @param fast logical. Should the Fast SVD be used? See details.
#' @param tol tolerance value below which the singular values of \code{x} are
#' considered zero.
#' @return a list containing the eigendecomposition for each class. If
#' \code{pool = TRUE}, then a single list is returned.
#' @examples
#' cov_eigen(x = iris[, -5], y = iris[, 5])
#' cov_eigen(x = iris[, -5], y = iris[, 5], pool = TRUE)
#' cov_eigen(x = iris[, -5], y = iris[, 5], pool = TRUE, fast = TRUE)
#'
#' # Generates a data set having fewer observations than features.
#' # We apply the Fast SVD to compute the eigendecomposition corresponding to the
#' # nonzero eigenvalues of the covariance matrices.
#' set.seed(42)
#' n <- 5
#' p <- 20
#' num_classes <- 3
#' x <- lapply(seq_len(num_classes), function(k) {
#'   replicate(p, rnorm(n, mean = k))
#' })
#' x <- do.call(rbind, x)
#' y <- gl(num_classes, n)
#' cov_eigen(x = x, y = y, fast = TRUE)
#' cov_eigen(x = x, y = y, pool = TRUE, fast = TRUE)
cov_eigen <- function(x, y, pool = FALSE, fast = FALSE, tol = 1e-6) {
  x <- as.matrix(x)
  y <- as.factor(y)

  if (fast) {
    # First, we center each class.
    x_centered <- center_data(x = x, y = y)
    if (pool) {
      n <- nrow(x_centered)
      svd_out <- fast.svd(x_centered)
      eigenvals <- svd_out$d^2 / n
      # We retain only the 'q' eigenvalues greater than 'tol'. Effectively, 'q'
      # is our estimate for the rank of the covariance matrix.
      q <- sum(eigenvals > tol)
      seq_q <- seq_len(q)
      eigen_out <- list(vectors = svd_out$v[, seq_q], values = eigenvals[seq_q])
    }
    else {
      eigen_out <- tapply(seq_along(y), y, function(i) {
        n <- length(i)
        svd_out <- fast.svd(x_centered[i, ])
        eigenvals <- svd_out$d^2 / n
        # We retain only the 'q' eigenvalues greater than 'tol'. Effectively, 'q'
        # is our estimate for the rank of the covariance matrix.
        q <- sum(eigenvals > tol)
        seq_q <- seq_len(q)
        list(vectors = svd_out$v[, seq_q], values = eigenvals[seq_q])
      })
    }
  } else {
    if (pool) {
      eigen_out <- eigen(cov_pool(x = x, y = y), symmetric = TRUE)
    } else {
      eigen_out <- lapply(cov_list(x = x, y = y), eigen, symmetric = TRUE)
    }
  }
  eigen_out
}

#' Computes the observation weights for each class for the HDRDA classifier
#'
#' This function calculates the weight for each observation in the data matrix
#' \code{x} in order to calculate the covariance matrices employed in the HDRDA
#' classifier, implemented in \code{\link{hdrda}}.
#'
#' @param x matrix containing the training data. The rows are the sample
#' observations, and the columns are the features.
#' @param y vector of class labels for each training observation
#' @param lambda the RDA pooling parameter. Must be between 0 and 1, inclusively.
#' @return list containing the observations for each class given in \code{y}
#'
#' @references Ramey, J. A., Stein, C. K., and Young, D. M. (2013),
#' "High-Dimensional Regularized Discriminant Analysis."
rda_weights <- function(x, y, lambda = 1) {
  x <- as.matrix(x)
  y <- as.factor(y)
  N <- length(y)
  tapply(seq_along(y), y, function(i) {
    n_k <- length(i)
    weights_k <- lambda / N + (1 - lambda) / n_k
    replace(rep(lambda / N, N), i, weights_k)
  })
}

#' Calculates the RDA covariance-matrix estimators for each class
#'
#' For the classes given in the vector \code{y}, this function calculates the
#' class covariance-matrix estimators employed in the HDRDA classifier,
#' implemented in \code{\link{hdrda}}.
#'
#' @param x matrix containing the training data. The rows are the sample
#' observations, and the columns are the features.
#' @param y vector of class labels for each training observation
#' @param lambda the RDA pooling parameter. Must be between 0 and 1, inclusively.
#' @return list containing the RDA covariance-matrix estimators for each class
#' given in \code{y}
#'
#' @references Ramey, J. A., Stein, C. K., and Young, D. M. (2013),
#' "High-Dimensional Regularized Discriminant Analysis."
rda_cov <- function(x, y, lambda = 1) {
  x_centered <- center_data(x = x, y = y)
  weights <- rda_weights(x = x, y = y, lambda = lambda)
  lapply(weights, function(weights_k) {
    x_scaled <- sqrt(weights_k) * x_centered
    crossprod(x_scaled)
  })
}

#' Computes a shrunken version of the maximum likelihood estimator for the
#' sample covariance matrix under the assumption of multivariate normality.
#'
#' For a sample matrix, \code{x}, we compute the sample covariance matrix as the
#' maximum likelihood estimator (MLE) of the population covariance matrix and
#' shrink it towards its diagonal.
#'
#' Let \eqn{\widehat{\Sigma}} be the MLE of the covariance matrix \eqn{\Sigma}.
#' Then, we shrink the MLE towards its diagonal by computing
#' \deqn{\widehat{\Sigma}(\gamma) = \gamma \widehat{\Sigma} + (1 - \gamma)
#' \widehat{\Sigma} \circ I_p,} where \eqn{\circ} denotes the Hadamard product
#' and \eqn{\gamma \in [0,1]}.
#'
#' For \eqn{\gamma < 1}, the resulting shrunken covariance matrix estimator is
#' positive definite, and for \eqn{\gamma = 1}, we simply have the MLE, which can
#' potentially be positive semidefinite (singular).
#'
#' The estimator given here is based on Section 18.3.1 of the Hastie et
#' al. (2008) text.
#'
#' @export
#' @param x data matrix with \code{n} observations and \code{p} feature vectors
#' @param gamma the shrinkage parameter. Must be between 0 and 1, inclusively.
#' By default, the shrinkage parameter is 1, which simply yields the MLE.
#' @return shrunken sample covariance matrix of size \eqn{p \times p}
#'
#' @references Hastie, T., Tibshirani, R., and Friedman, J. (2008), "The
#' Elements of Statistical Learning: Data Mining, Inference, and Prediction,"
#' 2nd edition. \url{http://www-stat.stanford.edu/~tibs/ElemStatLearn/}
cov_shrink_diag <- function(x, gamma = 1) {
  if (gamma < 0 || gamma > 1) {
    stop("The value of 'gamma' must be between 0 and 1, inclusively.")
  }
  cov_x <- cov_mle(x)
  gamma * cov_x + (1 - gamma) * diag(diag(cov_x))
}
