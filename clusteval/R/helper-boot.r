#' Creates a list of indices for a stratified nonparametric bootstrap.
#'
#' This function creates a list of indices for a stratified nonparametric
#' bootstrap. Corresponding to our Cluster Omission Stability statistic
#' implemented in \code{\link{clustomit}}, we omit each group in turn and perform
#' a stratified bootstrap without the group. We denote the number of groups
#' as \code{num_clusters}, which is equal to \code{nlevels(factor(y))}.
#' Specifically, suppose that we omit the \eqn{k}th group. That is, we ignore all
#' of the observations corresponding to group \eqn{k}. Then, we sample with
#' replacement from each of the remaining groups (i.e., every group except for
#' group \eqn{k}), yielding a set of bootstrap indices.
#'
#' The returned list contains \eqn{K \times num_reps} elements.
#'
#' @export
#' @param y a vector that denotes the grouping of each observation. It must be
#' coercible with \code{as.factor}.
#' @param num_reps the number of bootstrap replications to use for each group
#' @return named list containing indices for each bootstrap replication
#' @examples
#' set.seed(42)
#' # We use 4 clusters, each with up to 10 observations. The sample sizes are
#' # randomly chosen.
#' num_clusters <- 4
#' sample_sizes <- sample(10, num_clusters, replace = TRUE)
#'
#' # Create the cluster labels, y.
#' y <- unlist(sapply(seq_len(num_clusters), function(k) {
#'  rep(k, sample_sizes[k])
#' }))
#'
#' # Use 20 reps per group.
#' boot_stratified_omit(y, num_reps = 20)
#'
#' # Use the default number of reps per group.
#' boot_stratified_omit(y)
boot_stratified_omit <- function(y, num_reps = 50) {
  y <- as.factor(y)

  # A sequence that indexes each bootstrap replication
  seq_reps <- seq_len(num_reps)

  # Creates a list with each named element containing the indices of its cluster
  # members
  y_index <- split(seq_along(y), y)

  # We create a list, where each element contains the indices (rows) to apply the
  # clustering algorithm to examine its stability. This is effectively our
  # approach to stratified sampling. Once we have computed a similarity score for
  # each bootstrap rep, we can determine the corresponding omitted cluster with
  # the 'unsplit' function.
  boot_indices <- lapply(seq_len(nlevels(y)), function(omit_k) {
    lapply(seq_reps, function(b) {
      unlist(lapply(y_index[-omit_k], sample, replace = TRUE), use.names = FALSE)
    })
  })
  boot_indices <- unlist(boot_indices, recursive = FALSE)
  names(boot_indices) <- paste("Rep", seq_along(boot_indices), sep = "")
  boot_indices
}
