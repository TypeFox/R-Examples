#' @title Select Tuning Parameter for Thresholding Covariance Matrix by CV
#'
#' @description
#' Apply K-fold cross-validation for selecting tuning parameters for
#' thresholding covariance matrix using grid search strategy
#'
#' @details
#' For cross-validation, this function split the sample randomly into
#' two pieces of size n1 = n-n/log(n) and n2 = n/log(n), and repeat this k times
#'
#' @param matrix a N*p matrix, N indicates sample size
#'   and p indicates the dimension
#' @param method thresholding method, "hard" or "soft"
#' @param thresh.len the number of thresholding values tested in
#'   cross-validation, the thresholding values will be a sequence of
#'   \code{thresh.len} equally spaced values from minimum threshold constant
#'   to largest covariance in sample covariance matrix
#' @param n.cv times that cross-validation repeated, the default number is 10
#' @param norm the norms used to measure the cross-validation errors,
#'   which can be the Frobenius norm "F" or the operator norm "O"
#' @param seed random seed, the default value is 142857
#' @return An object of class "CovCv" containing the cross-validation's result
#'   for covariance matrix regularization, including:
#'   \item{regularization}{regularization method, which is "Hard Thresholding"
#'     or "Soft Thresholding"}
#'   \item{parameter.opt}{selected optimal parameter by cross-validation}
#'   \item{cv.error}{the corresponding cross-validation errors}
#'   \item{n.cv}{times that cross-validation repeated}
#'   \item{norm}{the norm used to measure the cross-validation error}
#'   \item{seed}{random seed}
#'   \item{threshold.grid}{thresholding values tested in cross-validation}
#' @references "High-Dimensional Covariance Estimation" by Mohsen Pourahmadi
#' @examples
#' data(m.excess.c10sp9003)
#' retcov.cv <- threshold.cv(m.excess.c10sp9003, method = "hard",
#'                           thresh.len = 20, n.cv = 10, norm = "F", seed = 142857)
#' summary(retcov.cv)
#' plot(retcov.cv)
#' # Low dimension
#' @importFrom stats cov
#' @export

threshold.cv <- function(matrix, method = "hard", thresh.len = 20, n.cv = 10, norm = "F", seed = 142857) {
  if ((method %in% c("hard","soft")) == FALSE) {
    stop("This function only support two thresholding methods: hard and soft")
  }
  if ((norm %in% c("F","O","f","o")) == FALSE) {
    stop("This function only support two norm: Frobenius and operator")
  }
  # Threshold Loss Function
  thresh.loss <- function(mat1, mat2, threshold, method, norm) {
    cov.mat1 <- cov(mat1)
    cov.mat2 <- cov(mat2)
    if (method == "hard") {
      mat.diff <- hard.thresholding(cov.mat1, threshold) - cov.mat2
    } else if (method == "soft") {
      mat.diff <- soft.thresholding(cov.mat1, threshold) - cov.mat2
    }
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
  # Threshold Value Lower and Upper Bound
  cov.SAM <- cov(matrix)
  min.t <- threshold.min(cov(matrix), method = method)
  max.t <- max(abs(cov.SAM[cov.SAM != diag(diag(cov.SAM))]))
  # Threshold Values for Cross-Validation
  thresh.grid <- seq(from = min.t, to = max.t, length.out = thresh.len)
  # Cross-Validation
  loss.mat <- matrix(0, nrow = n.cv, ncol = thresh.len)
  for (i in 1:n.cv) {
    set.seed(seed + i)
    index <- sample(1:N, size = n1, replace = FALSE)
    mat1 <- matrix[index,]
    mat2 <- matrix[-index,]
    loss.mat[i,] <- sapply(thresh.grid, FUN = thresh.loss, mat1 = mat1, mat2 = mat2, method = method, norm = norm)
  }
  loss.vec <- colMeans(loss.mat)
  threshold.opt <- thresh.grid[which.min(loss.vec)]
  if (method == "hard") {
    regularization <- "Hard Thresholding"
  }else if(method == "soft") {
    regularization <- "Soft Thresholding"
  }
  result <- list(regularization = regularization,
                 parameter.opt = threshold.opt,
                 cv.error = loss.vec,
                 n.cv = n.cv, norm = norm, seed = seed,
                 threshold.grid = thresh.grid
  )
  class(result) <- "CovCv"
  return(result)
}
