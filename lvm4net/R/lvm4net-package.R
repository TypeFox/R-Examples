#' Latent Variable Models for Networks
#'
#' \code{lvm4net} provides a range of tools for latent variable models for
#' network data. Most of the models are implemented using a fast
#' variational inference approach.
#' Latent space models for binary networks: the function \code{\link{lsm}} implements the latent space model (LSM) introduced by  Hoff et al. (2002) using a variational inference and squared Euclidian distance; the function
#' \code{\link{lsjm}} implements latent space joint model (LSJM) for multiplex networks introduced by
#' Gollini and Murphy (2014).
#' These models assume that each node of a network has a latent position
#' in a latent space: the closer two nodes are in the latent space, the more likely
#' they are connected.
#' Functions for binary bipartite networks will be added soon.
#' @references Gollini, I., and Murphy, T. B. (2014), "Joint Modelling of Multiple Network Views", Journal of Computational and Graphical Statistics \url{http://arxiv.org/abs/1301.3759}.
#' @references Hoff, P., Raftery, A., and Handcock, M. (2002), "Latent Space Approaches to Social Network Analysis", Journal of the American Statistical Association, 97, 1090--1098.
#'
#' @name lvm4net-package
#' @aliases lvm4net
#' @import MASS
#' @import ergm
#' @import network
#' @import ellipse
#' @importFrom igraph layout.fruchterman.reingold graph.adjacency
#' @docType package
NULL
