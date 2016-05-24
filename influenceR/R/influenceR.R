
#' influenceR: Software tools to quantify structural importance of nodes in a network.
#'
#' The influenceR package includes functions to quantify the structural
#' importance of nodes in a network. Algorithms include Betweenness Centrality,
#' Bridging, Constraint Index, Effective Network Size, and Key Players.
#' Currently, algorithms are only guaranteed to work on undirected graphs; work
#' on directed graphs is in progress. These functions run on graph objects from
#' the igraph package.
#' 
#' In addition to igraph, this package makes use of the SNAP framework for a
#' high-performance graph data structure and an OpenMP-parallelized
#' implementation of Betweenness Centrality. See
#' \url{http://snap-graph.sourceforge.net}
#' 
#' @section Funding:
#' Development of this software package was supported by NIH grant R01 DA033875.
#' 
#' @section References:
#' The website and source code is located at \url{http://github.com/rcc-uchicago/influenceR}.
#'
#' @docType package
#' @name influenceR
#' @useDynLib influenceR
NULL
#> NULL