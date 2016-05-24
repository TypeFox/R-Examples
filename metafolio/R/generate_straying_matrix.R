#' Generate a matrix of straying proportions within a metapopulation
#'
#' Generate a matrix of straying proportions within a metapopulation.
#' Based on Eq. 2 in Cooper and Mangel (1999).
#'
#' @param n_pop Number of subpopulations.
#' @param stray_fraction Fraction of individuals that stray from a given
#'   subpopulation.
#' @param stray_decay_rate Exponential rate that straying decays with distance
#'   between subpopulations.
#' @export
#' @references
#' Cooper, A.B. and Mangel, M. 1999. The dangers of ignoring
#' metapopulation structure for the conservation of salmonids. Fish.
#' Bull. 97(2): 213-226.
#' @examples
#' x <- generate_straying_matrix(10, 0.01, 0.3)
#' image(x, col = rev(heat.colors(12)))
generate_straying_matrix <- function(n_pop, stray_fraction,
  stray_decay_rate){

  stray_mat <- matrix(ncol = n_pop, nrow = n_pop)

  # distribute the fraction of fish among the other possible
  # populations:
  for(i in 1:n_pop) {
    for(j in 1:n_pop) {
      # Cooper and Mangel 1999 eq (2):
      stray_mat[i,j] <- exp(-stray_decay_rate * abs(i - j))
    }
  }
  diag(stray_mat) <- 0 # no self-straying
  # so no salmon or lost or gained:
  normalization <- rowSums(stray_mat)
  # to match the desired overall straying rate by stream:
  stray_mat <- (stray_mat / normalization) * stray_fraction
  return(stray_mat)
}
