#' Snowball Sampling with Multiple Inclusion Around a Single Seed
#'
#' This function performs snowball sampling with multiple inclusions (LSMI) around a
#' single seed.
#'
#' @note \code{\link{LSMI}.}
#'
#' @param seed0 \code{num}. Id of a seed to be sampled around.
#' @references Thompson, M. E., Ramirez Ramirez, L. L., Lyubchich, V. and
#'  Gel, Y. R. (2015), Using the bootstrap for statistical inference
#'  on random graphs. Can J Statistics. doi: 10.1002/cjs.11271
#' @inheritParams Oempdegreedistrib
#' @return a list containing:
#'    \item{seed}{seed0 \code{num}. Id of a seed to be sampled around.}
#'    \item{sampleN}{A vector of numeric ids of the nodes from
#'          LSMI along with the original seed. This vector may have
#'          duplicates, since the algorithm allows for multiple inclusions.}
#'    \item{unode}{A vector containing the unique values in \code{$sampleN}.}
#'    \item{nodes.waves}{A list of length \code{n.neigh} containing vectors where
#'          each vector reports numeric ids of nodes sampled in a particular wave.}
#' @export
#' @examples
#' net <- artificial_networks[[1]]
#' a <- sample_about_one_seed(net = net, seed0 = 1, n.neigh = 2)
sample_about_one_seed <- function(net, seed0, n.neigh = 1) {

      sampleN <- nodes <- seed0
      nodes.waves <- as.list(rep(0, n.neigh))
      effEdges <- net$edges
      more <- TRUE
      nn <- n.neigh
      new.nodes <- 0

      wave <- 1
      while (wave <= n.neigh & more) {
            a <- is.element(effEdges, nodes)
            # ^ 'nodes' will be accumulating all included nodes (non repeated)
            if (any(a))
            {
                  eedges <- which(matrix(a, dim(effEdges)[1], 2), arr.ind = TRUE)
                  #^ now it is the row number and column where they are in the edges matrix
                  nodes.waves[[wave]] <- arr.nodes <- sort(effEdges[cbind(eedges[, 1], sapply(eedges[, 2], FUN = switch, 2,
                                                                                              1))])
                  # ^ the nodes we arrived to (duplicity is allowed)
                  # I need this specially to know which nodes were the last added:
                  if (!anyDuplicated(eedges[, 1])) {
                        new.nodes <- arr.nodes
                        #^ all the nodes we arrive to weren't already included in 'nodes'
                  } else {
                        new.nodes <- setdiff(arr.nodes, nodes)
                  }  #Then, already included nodes are not considered new because of inclusion of edge connecting them
                  ### subEdges<-effEdges[unique(a),] #the subset of edges. The repeated just have to be included once. maybe we are arriving
                  ### to the nodes more than once (due to small cycles) or we can get again to already included nodes (due to larger
                  ### cycles).  We want to include them as may times as they are neighbours of already included nodes. That is why I
                  ### consider arr.nodes.if a originally seed vertex is included more than once, it is because it was selected also by
                  ### following one edge and then it also has the category of non seed.

                  sampleN <- sort(c(sampleN, arr.nodes))
                  nodes <- unique(sampleN)
                  if (nn > 1)
                        effEdges <- effEdges[-unique(eedges[, 1]), ]  #I remove the 'used edges' to facilitate following searches within while,and assure
                  # we do not 'arrive' to a node more times than edges it has.
                  if (length(effEdges) > 0) {
                        if (is.vector(effEdges))
                              effEdges <- t(effEdges)
                        #I have to this very often in R to make sure effEdges is a matrix and not a vector
                  } else {
                        more <- FALSE
                  }  #when it reduces to become a matrix with one row.
            }  #end if(any(a))
            wave <- wave + 1
      }  #end while
      list(seed = seed0, sampleN = sampleN, unodes = nodes, nodes.waves = nodes.waves)
}
