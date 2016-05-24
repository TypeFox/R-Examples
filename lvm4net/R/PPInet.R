#' PPI genetic and physical interactions data
#' 
#' The dataset contains two undirected networks formed by genetic and
#' physical protein-protein interactions (PPI) between 67 Saccharomyces 
#' cerevisiae proteins.
#' The genetic interactions network is formed of 294 links, and the physical
#' interactions network is formed of 190 links.
#' The data were downloaded from the Biological General Repository for
#' Interaction Datasets (BioGRID) database \url{http://thebiogrid.org/}
#' 
#' \itemize{
#' \item \code{PPIgen} Binary adjacency matrix containing genetic interactions between 67 proteins.
#' \item \code{PPIphy} Binary adjacency matrix containing physical interactions between 67 proteins.
#' }
#' 
#' @references Gollini, I., and Murphy, T. B. (2014), "Joint Modelling of Multiple Network Views", Journal of Computational and Graphical Statistics \url{http://arxiv.org/abs/1301.3759}.
#' @docType data
#' @keywords datasets
#' @seealso \code{\link{PPIgen}}, \code{\link{PPIphy}}
#' @format Two binary adjacency matrices
#' @source \url{http://thebiogrid.org/}
#' @name PPInet
NULL
