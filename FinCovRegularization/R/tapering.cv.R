#' @title Select Tuning Parameter for Tapering Covariance Matrix by CV
#'
#' @description
#' Apply K-fold cross-validation for selecting tuning parameters for
#' tapering covariance matrix using grid search strategy
#'
#' @details
#' For cross-validation, this function split the sample randomly into
#' two pieces of size n1 = n-n/log(n) and n2 = n/log(n), and repeat this k times
#'
#' @param matrix a N*p matrix, N indicates sample size
#'   and p indicates the dimension
#' @param h the ratio between taper l_h and parameter l
#' @param n.cv times that cross-validation repeated, the default number is 10
#' @param norm the norms used to measure the cross-validation errors,
#'   which can be the Frobenius norm "F" or the operator norm "O"
#' @param seed random seed, the default value is 142857
#' @return An object of class "CovCv" containing the cross-validation's result
#'   for covariance matrix regularization, including:
#'   \item{regularization}{regularization method, which is "Tapering"}
#'   \item{parameter.opt}{selected optimal parameter by cross-validation}
#'   \item{cv.error}{the corresponding cross-validation errors}
#'   \item{n.cv}{times that cross-validation repeated}
#'   \item{norm}{the norm used to measure the cross-validation error}
#'   \item{seed}{random seed}
#' @references "High-Dimensional Covariance Estimation" by Mohsen Pourahmadi
#' @examples
#' data(m.excess.c10sp9003)
#' retcov.cv <- tapering.cv(m.excess.c10sp9003, n.cv = 10,
#'                          norm = "F", seed = 142857)
#' summary(retcov.cv)
#' plot(retcov.cv)
#' # Low dimension
#' @importFrom stats cov
#' @export

tapering.cv <- function(matrix, h = 1/2, n.cv = 10, norm = "F", seed = 142857) {
  if ((norm %in% c("F","O","f","o")) == FALSE) {
    stop("This function only support two norm: Frobenius and operator")
  }
  # Tapering Loss Function
  tapering.loss <- function(mat1, mat2, l, h, norm) {
    cov.mat1 <- cov(mat1)
    cov.mat2 <- cov(mat2)
    mat.diff <- tapering(cov.mat1, l, h) - cov.mat2
    if (norm %in% c("F","f")) {
      loss <- F.norm2(mat.diff)
    } else if(norm %in% c("O","o")) {
      loss <- O.norm2(mat.diff)
    }
    return(loss)
  }
  # Numbers for Splitting
  N <- dim(matrix)[1]
  n1 <- ceiling(N*(1 - 1/log(N)))
  n2 <- floor(N/log(N))
  # Tapering Values for Cross-Validation
  l.grid <- 0:(ncol(matrix) - 1)
  # Cross-Validation
  loss.mat <- matrix(0, nrow = n.cv, ncol = length(l.grid))
  for (i in 1:n.cv) {
    set.seed(seed + i)
    index <- sample(1:N, size = n1, replace = FALSE)
    mat1 <- matrix[index,]
    mat2 <- matrix[-index,]
    loss.mat[i,] <- sapply(l.grid, FUN = tapering.loss, mat1 = mat1, mat2 = mat2, h = h, norm = norm)
  }
  loss.vec <- colMeans(loss.mat)
  l.opt <- l.grid[which.min(loss.vec)]
  result <- list(regularization = "Tapering",
                 parameter.opt = l.opt,
                 cv.error = loss.vec,
                 n.cv = n.cv, norm = norm, seed = seed,
                 h = h)
  class(result) <- "CovCv"
  return(result)
}
