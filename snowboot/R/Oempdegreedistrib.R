#' Obtaining an Empirical Network Degree Distribution from Labeled Snowball
#' Sampling with Multiple Inclusion (LSMI).
#'
#' Oempdegreedistrib is used to obtain the empirical network degree distribution
#' from labeled snowball sampling with multiple inclusion (LSMI).
#' @references Thompson, M. E., Ramirez Ramirez, L. L., Lyubchich, V. and
#'    Gel, Y. R. (2015), Using the bootstrap for statistical inference
#'    on random graphs. Can J Statistics. doi: 10.1002/cjs.11271
#' @param net A network object that is list containing:
#'  \describe{
#'    \item{edges}{The edgelist of the network. A two column
#'      \code{matrix} where each row is an edge.}
#'    \item{degree}{The degree sequence of the network, which is
#'      an \code{integer} vector of length n.}
#'    \item{n}{The network order.}
#'  }
#'    The object can be created by \code{\link{local.network.MR.new5}} or
#'    it can be imported.
#' @param n.seeds A number of seeds in the snowball sample.
#'    It must be a positive integer.
#' @param n.neigh A number of waves to be sampled around each seed in LSMI.
#'    For example, n.neigh = 0 corresponds to seeds only, and n.neigh = 1
#'    corresponds to sampling seeds and their first neighbors).
#'    Note that the algorithm allows for multiple inclusions.
#' @param num.sam A number for the LSMI repetitions. Default value is one.
#' @param seeds A matrix of dimension \code{num.sam} x \code{n.seeds} containing the
#'    numeric ids of the seeds to initiate sampling. Each row of the matrix
#'    corresponds to one LSMI sample. Note that this is an optional parameter.
#'    WARNING: As of now, this feature is only supported when
#'    parameter \code{n.neigh} is greater than zero.
#' @return A list consisting of
#'    \item{samples}{A list of length \code{num.sam} where each element
#'          is a list containing three tables:
#'          the frequency of degrees sampled from seeds,
#'          non-seeds including duplicated nodes,
#'          and non-seeds without duplications.}
#'    \item{values}{A list of length \code{num.sam} where each element is a
#'          vector containing the unique degree values sampled in each LSMI.}
#'    \item{Oemp}{A list of length \code{num.sam} where each element
#'          is a list containing two tables based on different methods to
#'          estimating the empirical distribution from the network sample
#'          (One list per LSMI).}
#'    \item{num.sam}{A number for the LSMI repetitions.}
#'    \item{val.seeds}{A list of length \code{num.sam} where each element
#'          is a vector of unique degree values sampled solely from seeds
#'          (One vector per LSMI).}
#'    \item{val.nonseed}{a list of length \code{num.sam} where each
#'          element is a vector of unique degree values sampled
#'          solely from non-seeds
#'          (One vector per LSMI).
#'          Note: This item is unreported when n.neigh equals zero.}
#'    \item{n.seeds}{the number of seeds in the snowball sample.}
#'    \item{n.neigh}{the number of waves carried out by the snowball sample. See
#'          input argument for details.}
#'    \item{p0.real}{proportion of nodes from the network with degree zero.
#'          Note: p0.real is unreported when n.neigh equals zero.}
#'    \item{p0.seeds}{a list of length \code{num.sam} where each
#'          element is the proportion of seeds with degree zero.
#'          (One element per LSMI).}
#'    \item{ekseed}{a list of length \code{num.sam} where each
#'          element is the sample mean of the seeds.
#'          Note that This is unreported when n.neigh equals zero.}
#'    \item{seeds1}{a matrix of dimension \code{num.sam} x \code{n.seeds} with
#'          the numeric seed ids. Each row corresponds to one LSMI.}
#'    \item{nodes_of_LSMI}{a list of length \code{num.sam} where each element is
#'          vector containing the numeric ids of the nodes sampled using LSMI
#'          (One element per LSMI). Note: nodes_of_LSMI is unreported when
#'          n.neigh equals zero.}
#' @export
#' @examples
#' net <- artificial_networks[[1]]
#' sam.out <- Oempdegreedistrib(net = net, n.seeds = 40, n.neigh = 1, num.sam = 1)

Oempdegreedistrib <- function(net, n.seeds, n.neigh, num.sam = 1, seeds = NULL) {
      if (n.neigh == 0) {
            # only information from the seeds
            res <- Oempdegreedistrib0(net, n.seeds, n.neigh, num.sam, seeds)
      } else {
            res <- OempdegreedistribK(net, n.seeds, n.neigh, num.sam, seeds)
      }
      res
}
