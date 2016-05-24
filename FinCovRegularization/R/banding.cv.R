#' @title Select Tuning Parameter for Banding Covariance Matrix by CV
#'
#' @description
#' Apply K-fold cross-validation for selecting tuning parameters for
#' banding covariance matrix using grid search strategy
#'
#' @details
#' For cross-validation, this function split the sample randomly into
#' two pieces of size n1 = n-n/log(n) and n2 = n/log(n), and repeat this k times
#'
#' @param matrix a N*p matrix, N indicates sample size
#'   and p indicates the dimension
#' @param n.cv times that cross-validation repeated, the default number is 10
#' @param norm the norms used to measure the cross-validation errors,
#'   which can be the Frobenius norm "F" or the operator norm "O"
#' @param seed random seed, the default value is 142857
#' @return An object of class "CovCv" containing the cross-validation's result
#'   for covariance matrix regularization, including:
#'   \item{regularization}{regularization method, which is "Banding"}
#'   \item{parameter.opt}{selected optimal parameter by cross-validation}
#'   \item{cv.error}{the corresponding cross-validation errors}
#'   \item{n.cv}{times that cross-validation repeated}
#'   \item{norm}{the norm used to measure the cross-validation error}
#'   \item{seed}{random seed}
#' @references "High-Dimensional Covariance Estimation" by Mohsen Pourahmadi
#' @examples
#' data(m.excess.c10sp9003)
#' retcov.cv <- banding.cv(m.excess.c10sp9003, n.cv = 10,
#'                         norm = "F", seed = 142857)
#' summary(retcov.cv)
#' plot(retcov.cv)
#' # Low dimension
#' @importFrom stats cov
#' @export

banding.cv <- function(matrix, n.cv = 10, norm = "F", seed = 142857) {
  if ((norm %in% c("F","O","f","o")) == FALSE) {
    stop("This function only support two norm: Frobenius and operator")
  }
  # Banding Loss Function
  banding.loss <- function(mat1, mat2, k, norm) {
    cov.mat1 <- cov(mat1)
    cov.mat2 <- cov(mat2)
    mat.diff <- banding(cov.mat1, k) - cov.mat2
    if (norm %in% c("F","f")) {
      loss <- F.norm2(mat.diff)
    } else if (norm %in% c("O","o")) {
      loss <- O.norm2(mat.diff)
    }
    return(loss)
  }
  # Numbers for Splitting
  N <- dim(matrix)[1]
  n1 <- ceiling(N*(1 - 1/log(N)))
  n2 <- floor(N/log(N))
  # Banding Values for Cross-Validation
  k.grid <- 0:(ncol(matrix) - 1)
  # Cross-Validation
  loss.mat <- matrix(0, nrow = n.cv, ncol = length(k.grid))
  for (i in 1:n.cv) {
    set.seed(seed + i)
    index <- sample(1:N, size = n1, replace = FALSE)
    mat1 <- matrix[index,]
    mat2 <- matrix[-index,]
    loss.mat[i,] <- sapply(k.grid, FUN = banding.loss, mat1 = mat1, mat2 = mat2, norm = norm)
  }
  loss.vec <- colMeans(loss.mat)
  k.opt <- k.grid[which.min(loss.vec)]
  result <- list(regularization = "Banding",
                 parameter.opt = k.opt,
                 cv.error = loss.vec,
                 n.cv = n.cv, norm = norm, seed = seed
  )
  class(result) <- "CovCv"
  return(result)
}
