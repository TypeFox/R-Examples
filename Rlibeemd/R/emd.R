#' EMD decomposition
#' 
#' Decompose input data to Intrinsic Mode Functions (IMFs) with the
#' Empirical Mode Decomposition algorithm.
#'
#' This is a wrapper around \code{eemd} with \code{ensemble_size = 1} and \code{noise_strength = 0}.
#'
#' @export
#' @name emd
#' @param input Vector of length N. The input signal to decompose.
#' @param num_imfs Number of Intrinsic Mode Functions (IMFs) to compute. If num_imfs is set to zero, a value of
#'        num_imfs = emd_num_imfs(N) will be used, which corresponds to a maximal number of
#'        IMFs. Note that the final residual is also counted as an IMF in this
#'        respect, so you most likely want at least num_imfs=2.
#' @param S_number Integer. Use the S-number stopping criterion [1] for the EMD procedure with the given values of S.
#'        That is, iterate until the number of extrema and zero crossings in the
#'        signal differ at most by one, and stay the same for S consecutive
#'        iterations. Typical values are in the range 3--8. If \code{S_number} is
#'        zero, this stopping criterion is ignored. Default is 4.
#' @param num_siftings Use a maximum number of siftings as a stopping criterion. If
#'        \code{num_siftings} is zero, this stopping criterion is ignored. Default is 50.
#' @return Time series object of class \code{"mts"} where series corresponds to
#'        IMFs of the input signal, with the last series being the final residual.
#'  @references
#' \enumerate{
#'       \item{N. E. Huang, Z. Shen and S. R. Long, "A new view of nonlinear water
#'       waves: The Hilbert spectrum", Annual Review of Fluid Mechanics, Vol. 31
#'       (1999) 417--457}
#'       }
#' @seealso \code{\link{eemd}}, \code{\link{ceemdan}} 
emd <- function(input, num_imfs = 0, S_number = 4L, num_siftings = 50L) {
  if (!all(is.finite(input))) 
    stop("'input' must contain finite values only.")
  if (num_imfs < 0)
    stop("Argument 'num_imfs' must be non-negative integer.")
  if (S_number < 0)
    stop("Argument 'S_number' must be non-negative integer.")
  if (num_siftings < 0)
    stop("Argument 'num_siftings' must be non-negative integer.")
  
  output <- eemdR(input, num_imfs, ensemble_size = 1L, 
    noise_strength = 0L, S_number, num_siftings, 
    rng_seed = 0L, threads = 0L)
  if (inherits(input, "ts")) {
    tsp(output) <- tsp(input)
  } else tsp(output) <- c(1, nrow(output), 1)
  if (ncol(output) > 1) {
    class(output) <- c("mts", "ts", "matrix")
    colnames(output) <- c(paste("IMF", 1:(ncol(output) - 1)), "Residual")
  } else class(output) <- "ts"
  output
}