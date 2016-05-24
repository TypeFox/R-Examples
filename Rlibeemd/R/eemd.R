#' EEMD Decomposition
#' 
#' Decompose input data to Intrinsic Mode Functions (IMFs) with the Ensemble Empirical Mode 
#' Decomposition algorithm [1].
#' 
#' The size of the ensemble and the relative magnitude of the added noise are given by parameters 
#' \code{ensemble_size} and \code{noise_strength}, respectively.  The stopping criterion for the 
#' decomposition is given by either a S-number [2] or an absolute number of siftings. In the case 
#' that both are positive numbers, the sifting ends when either of the conditions is fulfilled.
#' 
#' @export
#' @name eemd
#' @param input Vector of length N. The input signal to decompose.
#' @param num_imfs Number of Intrinsic Mode Functions (IMFs) to compute. If num_imfs is set to zero,
#'   a value of num_imfs = emd_num_imfs(N) will be used, which corresponds to a maximal number of 
#'   IMFs. Note that the final residual is also counted as an IMF in this respect, so you most 
#'   likely want at least num_imfs=2.
#' @param ensemble_size Number of copies of the input signal to use as the ensemble.
#' @param noise_strength Standard deviation of the Gaussian random numbers used as additional noise.
#'   \bold{This value is relative} to the standard deviation of the input signal.
#' @param S_number Integer. Use the S-number stopping criterion for the EMD procedure with the given
#'   values of $S$. That is, iterate until the number of extrema and zero crossings in the signal 
#'   differ at most by one, and stay the same for S consecutive iterations. Typical values are in 
#'   the range 3--8. If \code{S_number} is zero, this stopping criterion is ignored. Default is 4.
#' @param num_siftings Use a maximum number of siftings as a stopping criterion. If 
#'   \code{num_siftings} is zero, this stopping criterion is ignored. Default is 50.
#' @param threads Non-negative integer defining the maximum number of parallel threads (via OpenMP's
#'   \code{omp_set_num_threads}. Default value 0 uses all available threads defined by OpenMP's 
#'   \code{omp_get_max_threads}.
#' @param rng_seed A seed for the GSL's Mersenne twister random number generator. A value of zero 
#'   (default) denotes an implementation-defined default value.
#' @return Time series object of class \code{"mts"} where series corresponds to IMFs of the input 
#'   signal, with the last series being the final residual.
#'   
#' @references \enumerate{ \item{Z. Wu and N. Huang, "Ensemble Empirical Mode Decomposition: A 
#'   Noise-Assisted Data Analysis Method", Advances in Adaptive Data Analysis, Vol. 1 (2009) 1--41} 
#'   \item{N. E. Huang, Z. Shen and S. R. Long, "A new view of nonlinear water waves: The Hilbert 
#'   spectrum", Annual Review of Fluid Mechanics, Vol. 31 (1999) 417--457} }
#' @seealso \code{\link{ceemdan}}
#' @examples
#' x <- seq(0, 2*pi, length.out = 500)
#' signal <- sin(4*x)
#' intermittent <- 0.1 * sin(80 * x)
#' y <- signal * (1 + ifelse(signal > 0.7, intermittent, 0))
#' 
#' plot(x = x,y = y,type = "l")
#' # Decompose with EEMD
#' imfs <- eemd(y, num_siftings = 10, ensemble_size = 50, threads = 1)
#'  
#' plot(imfs)
#' # High frequencies
#' ts.plot(rowSums(imfs[, 1:3]))
#' # Low frequencies
#' ts.plot(rowSums(imfs[, 4:ncol(imfs)]))
eemd <- function(input, num_imfs = 0, ensemble_size = 250L, 
  noise_strength = 0.2, S_number = 4L, num_siftings = 50L, 
  rng_seed = 0L, threads = 0L) {
  
  if (!all(is.finite(input))) 
    stop("'input' must contain finite values only.")
  if (num_imfs < 0)
    stop("Argument 'num_imfs' must be non-negative integer.")
  if (ensemble_size < 0)
    stop("Argument 'ensemble_size' must be non-negative integer.")
  if (noise_strength < 0)
    stop("Argument 'noise_strength' must be non-negative.")
  if (S_number < 0)
    stop("Argument 'S_number' must be non-negative integer.")
  if (num_siftings < 0)
    stop("Argument 'num_siftings' must be non-negative integer.")
  if (rng_seed < 0)
    stop("Argument 'rng_seed' must be non-negative integer.")
  if (threads < 0)
    stop("Argument 'threads' must be non-negative integer.")
  output <- eemdR(input, num_imfs, ensemble_size, 
    noise_strength, S_number, num_siftings, rng_seed, threads)
  if (inherits(input, "ts")) {
    tsp(output) <- tsp(input)
  } else tsp(output) <- c(1, nrow(output), 1)
  if (ncol(output) > 1) {
    class(output) <- c("mts", "ts", "matrix")
    colnames(output) <- c(paste("IMF", 1:(ncol(output) - 1)), "Residual")
  } else class(output) <- "ts"
  output
}
