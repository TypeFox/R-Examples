#' Convex Clustering via Splitting Methods
#' 
#' Clustering is a fundamental problem in science and engineering. Many classic methods such as \eqn{k}-means,
#' Gaussian mixture models, and hierarchical clustering, however, employ greedy algorithms which can be
#' entrapped in local minima, sometimes drastical suboptimal ones at that. Recently introduced convex relaxations
#' of \eqn{k}-means and hierarchical clustering shrink cluster centroids toward one another and ensure a unique global minimizer. 
#' This package provides two variable splitting methods
#' \itemize{
#' \item{Alternating Method of Multipliers (ADMM)}
#' \item{Alternating Minimization Algorithm (AMA)}
#' }
#' for solving this convex formulation of the clustering problem. We seek the centroids \eqn{u_i} that minimize
#' \deqn{
#' \frac{1}{2} \sum_i || x_i - u_i||_2^2 + \gamma \sum_l w_{l} ||u_{l1} - u_{l2} ||
#' }
#' Two penalty norms are currently supported: 1-norm and 2-norm.
#'
#' The two main functions are \code{\link[cvxclustr]{cvxclust_path_admm}} and \code{\link[cvxclustr]{cvxclust_path_ama}} which compute the cluster paths using
#' the ADMM and AMA methods respectively. The function \code{\link[cvxclustr]{cvxclust}} is a wrapper function that calls either 
#' \code{cvxclust_path_admm} or \code{cvxclust_path_ama} (the default) to perform the computation.
#' 
#' The functions \code{\link[cvxclustr]{kernel_weights}} and \code{\link[cvxclustr]{knn_weights}} can be used in sequence
#' to compute weights that can improve the quality of the clustering paths.
#' 
#' The typical usage consists of three steps:
#' \itemize{
#' \item Compute weights \code{w}.
#' \item Generate a geometrically increasing regularization parameter sequence. Unfortunately a closed form expression for the minimum amount of penalization to get complete coalescence is currently unknown.
#' \item Call \code{\link[cvxclustr]{cvxclust}} using the data \code{X}, weights \code{w}, and regularization parameter sequence \code{gamma}.
#' }
#'
#' Cluster assignments can also be retrieved from the solution to the convex clustering problem.
#' Both \code{cvxclust_path_admm} and \code{cvxclust_path_ama} output an object of class \code{cvxclustobject}.
#' A cluster assignment can be extracted in two steps:
#' \itemize{
#' \item Call \code{\link[cvxclustr]{create_adjacency}} to construct an adjacency matrix from the centroid differences variable \code{V}.
#' \item Call \code{\link[cvxclustr]{find_clusters}} to extract the connected components of the adjacency matrix.
#' }
#' 
#' @author Eric C. Chi, Kenneth Lange
#' @references Eric C. Chi and Kenneth Lange. Splitting Methods
#'   for Convex Clustering. Journal of Computational and Graphical Statistics, in press.
#'   \url{http://arxiv.org/abs/1304.0499}.
#' @name cvxclustr
#' @docType package
NULL
