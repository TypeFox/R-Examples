#' Computes the similarity between two clusterings of the same data set.
#'
#' For two clusterings of the same data set, this function calculates the 
#' similarity statistic specified of the clusterings from the comemberships of
#' the observations. Basically, the comembership is defined as the pairs of
#' observations that are clustered together.
#'
#' To calculate the similarity, we compute the 2x2 contingency table, consisting
#' of the following four cells:
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
#' Currently, we have implemented the following similarity statistics:
#' \itemize{
#'   \item Rand index
#'   \item Jaccard coefficient
#' }
#'
#' To compute the contingency table, we use the \code{\link{comembership_table}}
#' function.
#'
#' @export
#' @param labels1 a vector of \code{n} clustering labels
#' @param labels2 a vector of \code{n} clustering labels
#' @param similarity the similarity statistic to calculate
#' @param method the model under which the statistic was derived
#' @return the similarity between the two clusterings
#' @examples
#' # Notice that the number of comemberships is 'n choose 2'.
#' iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
#' iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
#' cluster_similarity(iris_kmeans, iris_hclust)
cluster_similarity <- function(labels1, labels2,
                               similarity = c("jaccard", "rand"),
                               method = "independence") {
	similarity <- match.arg(similarity)
  method <- match.arg(method)

  # Currently, we ignore the `method` argument and only use the similarity
  # statistics derived under an independence assumption.
  switch(similarity,
         jaccard = jaccard_indep(labels1, labels2),
         rand = rand_indep(labels1, labels2))
}

#' Computes the Jaccard similarity coefficient of two clusterings of the same
#' data set under the assumption that the two clusterings are independent.
#'
#' For two clusterings of the same data set, this function calculates the Jaccard
#' similarity coefficient of the clusterings from the comemberships of the
#' observations. Basically, the comembership is defined as the pairs of
#' observations that are clustered together.
#'
#' To calculate the Rand index, we compute the 2x2 contingency table, consisting
#' of the following four cells:
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
#' The Jaccard similarity coefficient is defined as:
#' \deqn{J = \frac{n_{11}}{n_{11} + n_{10} + n_{01}}}.
#'
#' In the special case that the Jaccard coefficient results in \eqn{0/0},
#' we define \eqn{J = 0}. For instance, this case can occur when both clusterings
#' consist of all singleton clusters.
#'
#' To compute the contingency table, we use the \code{\link{comembership_table}}
#' function.
#'
#' @param labels1 a vector of \code{n} clustering labels
#' @param labels2 a vector of \code{n} clustering labels
#' @return the Jaccard coefficient for the two sets of cluster labels (See
#' Details.)
#' @examples
#'\dontrun{
#' # We generate K = 3 labels for each of n = 10 observations and compute the
#' # Jaccard similarity coefficient between the two clusterings.
#' set.seed(42)
#' K <- 3
#' n <- 10
#' labels1 <- sample.int(K, n, replace = TRUE)
#' labels2 <- sample.int(K, n, replace = TRUE)
#' jaccard_indep(labels1, labels2)
#' 
#' # Here, we cluster the \code{\link{iris}} data set with the K-means and
#' # hierarchical algorithms using the true number of clusters, K = 3.
#' # Then, we compute the Jaccard similarity coefficient between the two
#' # clusterings.
#' iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
#' iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
#' jaccard_indep(iris_kmeans, iris_hclust)
#' }
jaccard_indep <- function(labels1, labels2) {
  com_table <- comembership_table(labels1, labels2)
  jaccard_out <- with(com_table, n_11 / (n_11 + n_10 + n_01))

  # In the case where 'labels1' and 'labels2' contain all singletons, the Jaccard
  # coefficient results in the expression 0 / 0, which yields a NaN value in R.
  # We define such cases as 0.
  if (is.nan(jaccard_out)) {
    warning("The two clusterings contain all singletons -- returning 0.")
    jaccard_out <- 0
  }
  jaccard_out
}

#' Computes the Rand similarity index of two clusterings of the same data set
#' under the assumption that the two clusterings are independent.
#'
#' For two clusterings of the same data set, this function calculates the Rand
#' similarity coefficient of the clusterings from the comemberships of the
#' observations. Basically, the comembership is defined as the pairs of
#' observations that are clustered together.
#'
#' To calculate the Rand index, we compute the 2x2 contingency table, consisting
#' of the following four cells:
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
#' The Rand similarity index is defined as:
#' \deqn{R = \frac{n_{11} + n_{00}}{n_{11} + n_{10} + n_{01} + n_{00}}}.
#'
#' To compute the contingency table, we use the \code{\link{comembership_table}}
#' function.
#'
#' @param labels1 a vector of \code{n} clustering labels
#' @param labels2 a vector of \code{n} clustering labels
#' @return the Rand index for the two sets of cluster labels
#' @examples
#'\dontrun{
#' # We generate K = 3 labels for each of n = 10 observations and compute the
#' # Rand similarity index between the two clusterings.
#' set.seed(42)
#' K <- 3
#' n <- 10
#' labels1 <- sample.int(K, n, replace = TRUE)
#' labels2 <- sample.int(K, n, replace = TRUE)
#' rand_indep(labels1, labels2)
#' 
#' # Here, we cluster the \code{\link{iris}} data set with the K-means and
#' # hierarchical algorithms using the true number of clusters, K = 3.
#' # Then, we compute the Rand similarity index between the two clusterings.
#' iris_kmeans <- kmeans(iris[, -5], centers = 3)$cluster
#' iris_hclust <- cutree(hclust(dist(iris[, -5])), k = 3)
#' rand_indep(iris_kmeans, iris_hclust)
#' }
rand_indep <- function(labels1, labels2) {
  com_table <- comembership_table(labels1, labels2)
  with(com_table, (n_11 + n_00) / (n_11 + n_10 + n_01 + n_00))
}
