#' 200 Simulated Networks of order 2000 with Polylogarithmic (0.1, 2)
#' Degree Distributions
#'
#' A list called "networks" containing 200 network objects of order 2000. These
#' networks were simulated using the polylogarithmic (aka Gutenberg-Richter law)
#' degree distribution (Newman et al., 2001; Newman, 2002) with parameters
#' \eqn{\delta = 0.1} and \eqn{\lambda = 2} as see in the following equations:
#' \deqn{f(k) = k^{-{\delta}}e^{-{k/{\lambda}}}/Li_{\delta}(e^{-{1/\lambda}})}{f(k)=k^-\delta exp(-k/\lambda )/Li[\delta](exp(-1/\lambda))}
#' \deqn{Li=\sum_{j=1}^{\infty} z^{-j}/{j^{\delta}}}{Li=\sum_{j=1}^{\infty} z^{-j}/{j^{\delta}}}
#' where \eqn{\lambda > 0}. Please see refence below for details (Thompson  et al, 2016).
#' @references Thompson, M. E., Ramirez Ramirez, L. L., Lyubchich, V. and
#' Gel, Y. R. (2015), Using the bootstrap for statistical inference
#' on random graphs. Can J Statistics. doi: 10.1002/cjs.11271
#' @format a list containing 200 network objects. Each network object is a list
#' with three elements:
#' \describe{
#'    \item{edges}{The edgelist of the network. A two column
#'      \code{matrix} where each row is an edge.}
#'    \item{degree}{The degree sequence of the network, which is
#'      an \code{integer} vector of length n.}
#'    \item{n}{The network order. The order for every network is 2000.}
#' }



"artificial_networks"
