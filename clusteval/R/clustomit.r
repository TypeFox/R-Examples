#' ClustOmit - Cluster Stability Evaluation via Cluster Omission
#'
#' We provide an implementation of the ClustOmit statistic, which is an approach
#' to evaluating the stability of a clustering determined by a clustering
#' algorithm. As discussed by Hennig (2007), arguably a stable clustering is one
#' in which a perturbation of the original data should yield a similar
#' clustering. However, if a perturbation of the data yields a large change in
#' the clustering, the original clustering is considered unstable. The ClustOmit
#' statistic provides an approach to detecting instability via a stratified,
#' nonparametric resampling scheme. We determine the stability of the clustering
#' via the similarity statistic specified (by default, the Jaccard coefficient).
#' 
#' To compute the ClustOmit statistic, we first cluster the data given in
#' \code{x} into \code{num_clusters} clusters with the clustering algorithm
#' specified in \code{cluster_method}. We then omit each cluster in turn
#' and all of the observations in that cluster. For the omitted cluster, we
#' resample from the remaining observations and cluster the resampled
#' observations into \code{num_clusters - 1} clusters again using the clustering
#' algorithm specified in \code{cluster_method}. Next, we compute the similarity
#' between the cluster labels of the original data set and the cluster labels of
#' the bootstrapped sample. We approximate the sampling distribution of the
#' ClustOmit statistic using a stratified, nonparametric bootstrapping scheme and
#' use the apparent variability in the approximated sampling distribution as a
#' diagnostic tool for further evaluation of the proposed clusters. By default,
#' we utilize the Jaccard similarity coefficient in the calculation of the
#' ClustOmit statistic to provide a clear interpretation of cluster assessment.
#' The technical details of the ClustOmit statistic can be found in our
#' forthcoming publication entitled "Cluster Stability Evaluation of Gene
#' Expression Data."
#'
#' The ClustOmit cluster stability statistic is based on the cluster omission
#' admissibility condition from Fisher and Van Ness (1971), who provide
#' decision-theoretic admissibility conditions that a reasonable clustering
#' algorithm should satisfy. The guidelines from Fisher and Van Ness (1971)
#' establish a systematic foundation that is often lacking in the evaluation of
#' clustering algorithms. The ClustOmit statistic is our proposed methodology to
#' evaluate the cluster omission admissibility condition from Fisher and
#' Van Ness (1971).
#'
#' We require a clustering algorithm function to be specified in the argument
#' \code{cluster_method}. The function given should accept at least two
#' arguments:
#' \describe{
#'   \item{x}{matrix of observations to cluster}
#'   \item{num_clusters}{the number of clusters to find}
#'   \item{...}{additional arguments that can be passed on}
#' }
#' Also, the function given should return only clustering labels for each
#' observation in the matrix \code{x}. The additional arguments specified in
#' \code{...} are useful if a wrapper function is used: see the example below for
#' an illustration.
#'
#' @export
#' @param x data matrix with \code{n} observations (rows) and \code{p} features
#' (columns)
#' @param num_clusters the number of clusters to find with the clustering
#' algorithm specified in \code{cluster_method}
#' @param cluster_method a character string or a function specifying the
#' clustering algorithm that will be used. The method specified is matched with
#' the \code{\link{match.fun}} function. The function given should return only
#' clustering labels for each observation in the matrix \code{x}.
#' @param similarity the similarity statistic that is used to compare the
#' original clustering (after a single cluster and its observations have been
#' omitted) to its resampled counterpart. Currently, we have implemented the
#' Jaccard and Rand similarity statistics and use the Jaccard statistic by
#' default.
#' @param weighted_mean logical value. Should the aggregate similarity score for
#' each bootstrap replication be weighted by the number of observations in each
#' of the observed clusters? By default, yes (i.e., \code{TRUE}).
#' @param num_reps the number of bootstrap replicates to draw for each omitted
#' cluster
#' @param num_cores the number of coures to use. If 1 core is specified, then
#' \code{\link{lapply}} is used without parallelization. See the \code{mc.cores}
#' argument in \code{\link{mclapply}} for more details.
#' @param ... additional arguments passed to the function specified in
#' \code{cluster_method}
#' @return object of class \code{clustomit}, which contains a named list with
#' elements
#' \describe{
#'   \item{boot_aggregate:}{vector of the aggregated similarity statistics for
#' each bootstrap replicate}
#'   \item{boot_similarity:}{list containing the bootstrapped similarity scores
#' for each cluster omitted}
#'   \item{obs_clusters:}{the clustering labels determined for the observations
#' in \code{x}}
#'   \item{num_clusters:}{the number of clusters found}
#'   \item{similarity:}{the similarity statistic used for comparison between the
#' original clustering and the resampled clusterings}
#' }
#' @references Fisher, L. and Van Ness, J. (1971), Admissible Clustering
#' Procedures, _Biometrika_, 58, 1, 91-104.
#' @references Hennic, C. (2007), Cluster-wise assessment of cluster stability,
#' _Computational Statistics and Data Analysis_, 52, 258-271.
#' \url{http://www.jstor.org/stable/2334320}
#' @examples
#' # First, we create a wrapper function for the K-means clustering algorithm
#' # that returns only the clustering labels for each observation (row) in
#' # \code{x}.
#' kmeans_wrapper <- function(x, num_clusters, num_starts = 10, ...) {
#'   kmeans(x = x, centers = num_clusters, nstart = num_starts, ...)$cluster
#' }
#'
#' # For this example, we generate five multivariate normal populations with the
#' # \code{sim_data} function.
#' x <- sim_data("normal", delta = 1.5, seed = 42)$x
#'
#' clustomit_out <- clustomit(x = x, num_clusters = 4,
#'                            cluster_method = "kmeans_wrapper", num_cores = 1)
#' clustomit_out2 <- clustomit(x = x, num_clusters = 5,
#'                             cluster_method = kmeans_wrapper, num_cores = 1)
clustomit <- function(x, num_clusters, cluster_method,
                      similarity = c("jaccard", "rand"), weighted_mean = TRUE,
                      num_reps = 50, num_cores = getOption("mc.cores", 2), ...) {
  x <- as.matrix(x)
  num_clusters <- as.integer(num_clusters)
  cluster_method <- match.fun(cluster_method)
  similarity <- match.arg(similarity)

  # The cluster labels for the observed (original) data matrix (i.e., x).  
  obs_clusters <- cluster_method(x = x, num_clusters = num_clusters, ...)
  cluster_sizes <- as.vector(table(obs_clusters))

  # Determines the indices for the bootstrap reps.
  boot_indices <- boot_stratified_omit(y = obs_clusters, num_reps = num_reps)

  # For each set of bootstrap indices, cluster the resampled data and
  # compute the similarity with the corresponding original clusters.
  boot_similarity <- mclapply(boot_indices, function(idx) {
    clusters_omit <- cluster_method(x = x[idx, ],
                                    num_clusters = num_clusters - 1, ...)

    cluster_similarity(obs_clusters[idx], clusters_omit,
                       similarity = similarity)
  }, mc.cores = num_cores)

  # Because 'mclapply' returns a list, we first simplify the list to an array and
  # then 'split' the similarity scores into a list by cluster.
  boot_similarity <- simplify2array(boot_similarity)
  boot_similarity <- split(boot_similarity, gl(num_clusters, num_reps))
  boot_similarity <- lapply(boot_similarity, as.vector)

  # Now, we compute the weighted average of the similarity scores for each
  # bootstrap replication. The weights correspond to the sample sizes of each
  # cluster.
  boot_similarity_matrix <- do.call(cbind, boot_similarity)
  boot_aggregate <- apply(boot_similarity_matrix, 1, weighted.mean,
                          w = cluster_sizes)

	obj <- list(
    boot_aggregate = as.vector(boot_aggregate),
		boot_similarity = boot_similarity,
    obs_clusters = obs_clusters,
		num_clusters = num_clusters,
		similarity = similarity
	)
	class(obj) <- "clustomit"
	obj
}
