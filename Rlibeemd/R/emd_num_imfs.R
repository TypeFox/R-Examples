#' Number of IMFs
#' 
#' Return the number of IMFs extracted from input data of length N, including
#' the final residual. This is just [log_2(N)] for N>3.
#' @export
#' @name nIMFs
#' @param N An integer defining the length of input data.
#' @return The number of IMFs which would be extracted from input data of length N, including
#' the final residual. 
emd_num_imfs <- function(N) {  
  if (!isTRUE(N > 0) || !isTRUE(abs(N - round(N)) < 100 * .Machine$double.eps))
    stop("N must be a positive integer.")
  emd_num_imfsR(N)
}