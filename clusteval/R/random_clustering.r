#' Randomly cluster a data set into K clusters.
#'
#' For each observation (row) in 'x', one of K labels is randomly generated.
#' By default, the probabilities of selecting each clustering label are equal,
#' but this can be altered by specifying 'prob', a vector of probabilities for
#' each cluster.
#'
#' Random clustering is often utilized as a baseline comparison clustering
#' against which other clustering algorithms are employed to identify structure
#' within the data. Typically, comparisons are made in terms of proposed
#' clustering assessment and evaluation methods as well as clustering similarity
#' measures. For the former, a specified clustering evaluation method is computed
#' for the considered clustering algorithms as well as random clustering. If the
#' clusters determined by a considered clustering algorithm do not differ
#' significantly from the random clustering, we might conclude that the found
#' clusters are no better than naively choosing clustering labels for each
#' observation at random. Likewise, a similarity measure can be computed to
#' compare the clusterings from each of a considered clustering algorithm and a
#' random clustering: if the clusterings are significantly similar, once again,
#' we might conclude the clusters found via the considered clustering algorithm
#' do not differ significantly from those found at random. In either case, the
#' clusters are unlikely to provide meaningful results on which the user can
#' better understand the inherent structure within the data.
#' 
#' @param x a matrix containing the data to cluster. The rows are the sample
#' observations, and the columns are the features.
#' @param K the number of clusters
#' @param prob a vector of probabilities to generate each cluster label. If
#' NULL, each cluster label has an equal chance of being selected.
#' @return a vector of clustering labels for each observation in 'x'.
random_clustering <- function(x, K, prob = NULL) {
  if (!is.null(prob)) {
    if (!is.numeric(prob)) {
      stop("The vector 'prob' must be 'numeric'.")
    }
    if (K != length(prob)) {
      stop("The length of 'prob' must equal 'K'.")
    }
    if (sum(prob) != 1) {
      stop("The sum of the probabilities must sum to 1.")
    }
    if (any(prob <= 0) || any(prob >= 1)) {
      stop("The cluster probabilties must be between 0 and 1.")
    }
  }

  sample(x = seq_len(K), size = nrow(x), replace = TRUE, prob = prob)
}
