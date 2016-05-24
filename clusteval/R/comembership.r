#' Calculates the comemberships of all pairs of a vector of clustering labels.
#'
#' For a set of clustering labels, this function computes the comembership of all
#' pairs of observations. Basically, two observations are said to be comembers if
#' they are clustered together.
#'
#' Tibshirani and Walther (2005) use the term 'co-membership', which we shorten
#' to 'comembership'. Some authors instead use the terms 'connectivity' or
#' 'co-occurrence'.
#' 
#' We use the \code{Rcpp} package to improve the runtime speed of this function.
#'
#' @export
#' @param labels a vector of \code{n} clustering labels
#' @return a vector of \code{choose(n, 2)} comembership bits
#' @references Tibshirani, R. and  Walther, G. (2005), Cluster Validation by
#' Prediction Strength, _Journal of Computational and Graphical Statistics_, 14,
#' 3, 511-528.
#' \url{http://amstat.tandfonline.com/doi/abs/10.1198/106186005X59243}.
#' @examples
#' # We generate K = 3 labels for each of n = 10 observations and compute the
#' # comembership for all 'n choose 2' pairs.
#' set.seed(42)
#' K <- 3
#' n <- 10
#' labels <- sample.int(K, n, replace = TRUE)
#' comembership_out <- comembership(labels)
#' comembership_out
#' 
#' # Notice that the number of comemberships is 'n choose 2'.
#' length(comembership_out) == choose(n, 2)
comembership <- function(labels) {
	.Call("rcpp_comembership", labels, PACKAGE = "clusteval")
}

#' Calculates the 2x2 contingency table of agreements and disagreements of
#' comemberships from two vectors of clustering labels.
#'
#' For two clusterings of the same data set, this function calculates the 2x2
#' contingency table of agreements and disagreements of the corresponding two
#' vectors of comemberships. Basically, the comembership is defined as the pairs
#' of observations that are clustered together.
#'
#' The contingency table calculated is typically utilized in the calculation of
#' a similarity statistic (e.g., Rand index, Jaccard index) between the two
#' clusterings. The 2x2 contingency table consists of the following four cells:
#' \describe{
#'   \item{n_11}{the number of observation pairs where both observations are
#' comembers in both clusterings}
#'   \item{n_10}{the number of observation pairs where the observations are
#' comembers in the first clustering but not the second}
#'   \item{n_01}{the number of observation pairs where the observations are
#' comembers in the second clustering but not the first}
#'   \item{n_00}{the number of observation pairs where neither pair are comembers
#' in either clustering}
#' }
#'
#' Tibshirani and Walther (2005) use the term 'co-membership', which we shorten
#' to 'comembership'. Some authors instead use the terms 'connectivity' or
#' 'co-occurrence'.
#' 
#' We use the \code{Rcpp} package to improve the runtime speed of this function.
#'
#' @export
#' @param labels1 a vector of \code{n} clustering labels
#' @param labels2 a vector of \code{n} clustering labels
#' @return named list containing the calculated contingency table:
#' \itemize{
#'   \item n_11
#'   \item n_10
#'   \item n_01
#'   \item n_00
#' }
#' @references Tibshirani, R. and  Walther, G. (2005). Cluster Validation by
#' Prediction Strength. Journal of Computational and Graphical Statistics, 14, 3,
#' 511-528. \url{http://amstat.tandfonline.com/doi/abs/10.1198/106186005X59243}.
#' @examples
#' # We generate K = 3 labels for each of n = 10 observations and compute the
#' # comembership for all 'n choose 2' pairs.
#' set.seed(42)
#' K <- 3
#' n <- 10
#' labels1 <- sample.int(K, n, replace = TRUE)
#' labels2 <- sample.int(K, n, replace = TRUE)
#' comembership_table(labels1, labels2)
#' 
#' # Here, we cluster the \code{\link{iris}} data set with the K-means and
#' # hierarchical algorithms using the true number of clusters, K = 3.
#' # Then, we compute the 2x2 contingency table agreements and disagreements of
#' #' the comemberships.
#' iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
#' iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
#' comembership_table(iris_kmeans, iris_hclust)
comembership_table <- function(labels1, labels2) {
  if (length(labels1) != length(labels2)) {
    stop("The two vectors of cluster labels must be of equal length.");
  }

	.Call("rcpp_comembership_table", labels1, labels2, PACKAGE = "clusteval")
}
