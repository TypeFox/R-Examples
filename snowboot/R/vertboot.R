#'Bootstrapping a network with vertex bootstrap
#'
#'This function bootstraps the the original network using a vertex bootstrap technique.
#'
#' @param m1 An adjacency matrix,the matrix represents the original network.
#' @param n.boot A positive integer number, the number of bootstrap replications.
#' @references Tom A.B.Snijders., Stephen P.Borgatti. (1999). Non-Parametric Standard Errors
#'and Tests for Network Statistics.
#'
#' @return A list of bootstrapped networks as adjacency matricies.
#' @export
#' @examples
#' graph_ex <- igraph::graph_from_edgelist(artificial_networks[[1]]$edges)
#' m1 <- igraph::as_adjacency_matrix(graph_ex)
#' m1 <- as.matrix(m1)
#' vertboot_out <- vertboot(m1,20)

vertboot <- function(m1, n.boot){
  res <- list()
  for (i in 1:n.boot) {
    blist <- sample(0:(dim(m1)[1]-1), replace = T)
    res <- c(res, list(vertboot_matrix_rcpp(m1,blist)))
  }
  res
}
