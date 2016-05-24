#' Extract the backbone of a weighted network using the disparity filter
#'
#' Given a weighted graph, \code{backbone} identifies the 'backbone structure'
#' of the graph, using the disparity filter algorithm by Serrano et al. (2009).
#' @param graph The input graph.
#' @param weights A numeric vector of edge weights, which defaults to
#' \code{E(graph)$weight}.
#' @param directed The directedness of the graph, which defaults to the result
#' of \code{\link[igraph]{is_directed}}.
#' @param alpha The significance level under which to preserve the edges, which
#' defaults to \code{0.05}.
#' @return An edge list corresponding to the 'backbone' of the graph, i.e. the
#' edges of the initial graph that were preserved by the null model that the
#' disparity filter algorithm implements.
#' @author Serrano et al. (2009); R implementation by Alessandro Bessi and
#' Francois Briatte
#' @references Serrano, M.A., Boguna, M. and Vespignani, A. (2009).
#' Extracting the multiscale backbone of complex weighted networks.
#' \emph{Proceedings of the National Academy of Sciences} 106, 6483-6488.
#' @examples
#' if (require(igraph)) {
#'
#'   # undirected network
#'   g <- sample_pa(n = 250, m = 5, directed = FALSE)
#'   E(g)$weight <- sample(1:25, ecount(g), replace = TRUE)
#'   backbone(g)
#'   
#'   # directed network
#'   g <- sample_pa(n = 250, m = 5, directed = TRUE)
#'   E(g)$weight <- sample(1:25, ecount(g), replace = TRUE)
#'   backbone(g)
#'
#' }
#' @aliases get.backbone
#' @importFrom igraph as_data_frame degree E ego is_directed is_igraph
#' @export
backbone <- function(graph, weights = igraph::E(graph)$weight,
                     directed = igraph::is_directed(graph), alpha = 0.05) {

  if (!igraph::is_igraph(graph)) {
    stop("Not a graph object")
  }

  stopifnot(!is.null(weights))

  if (!directed) {
    b = disparity_filter(graph, weights, "all", alpha)
  } else {
    b = rbind(
      disparity_filter(graph, weights, "in", alpha),
      disparity_filter(graph, weights, "out", alpha)
    )
  }
  return(unique(b[ order(b$from), ]))

}

#' @keywords internal
disparity_filter <- function(G, weights, mode = "all", alpha = 0.05) {

  d = igraph::degree(G, mode = mode)
  e = cbind(igraph::as_data_frame(G)[, 1:2 ], weight = weights, alpha = NA)
  if (mode == "all") {
    e = rbind(e, data.frame(from = e[, 2], to = e[, 1], e[, 3:4 ]))
  }

  for (u in which(d > 1)) {

    w = switch(substr(mode, 1, 1),
      a = which(e[, 1] == u | e[, 2] == u),
      i = which(e[, 2] == u),
      o = which(e[, 1] == u)
    )
    w = sum(e$weight[ w ]) / (1 + (mode == "all"))

    k = d[u]

    for (v in igraph::ego(G, 1, u, mode)[[1]][-1]) {

      ij = switch(substr(mode, 1, 1),
        a = which(e[, 1] == u & e[, 2] == v),
        i = which(e[, 1] == v & e[, 2] == u),
        o = which(e[, 1] == u & e[, 2] == v)
      )
      # cat(mode, "-", u, "->", v, ":", ij, "\n")

      # p_ij = e$weight[ ij ] / w
      # alpha_ij = integrate(function(x) { (1 - x) ^ (k - 2) }, 0, p_ij)
      # alpha_ij = 1 - (k - 1) * alpha_ij$value
      e$alpha[ ij ] = (1 - e$weight[ ij ] / w) ^ (k - 1)

    }

  }

  return(e[ !is.na(e$alpha) & e$alpha < alpha, 1:4 ])

}

# for compatibility with previous versions
get.backbone <- backbone
