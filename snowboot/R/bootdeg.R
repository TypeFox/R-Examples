#' Bootstrapping Empirical Degree Distribution
#'
#' This function delivers a bootstrap estimate of network degree distribution
#' based on a LSMI sample. Default is one bootstrap replication.
#'
#' @param sam.out A list that is the output of \code{\link{Oempdegreedistrib}}.
#' @param num.sam A vector of integers containing the numeric ids of the LSMI
#'    samples when \code{sam.out$num.sam} is greater than one. When \code{num.sam} is an
#'    integer, N, LSMI from 1 to N are taken from the input \code{sam.out}.
#' @param n.boot A positive integer number, the number of bootstrap replications.
#' @references Efron, B. (1979). Bootstrap methods: another look at the
#'  jackknife. The annals of Statistics, 1-26.
#' @references Thompson, M. E., Ramirez Ramirez, L. L., Lyubchich, V. and
#'  Gel, Y. R. (2015), Using the bootstrap for statistical inference
#'  on random graphs. Can J Statistics. doi: 10.1002/cjs.11271
#' @return A list consisting of:
#'    \item{values}{A list of length \code{num.sam} where each element is a
#'          vector containing the unique degree values sampled in each LSMI.}
#'    \item{empd}{A list of length three where each element contains a different
#'          estimate of degree distribution: empd.w.p0s - weighted bootstrap
#'          with a proportion of isolated nodes p0 being estimated by simple
#'          random sampling of bootstrapped seeds; empd.nw.p0sEkb - non-weighted
#'          bootstrap with a proportion of isolated nodes p0 being estimated by
#'          simple random sampling of bootstrapped seeds; empd.nw.p0sEks -ignore
#'          (see Thompson et al. for details)}
#'    \item{num.sam}{Numeric indices corresponding to LSMI samples used for bootstrap.}
#'    \item{n.boot}{The same object as input argument \code{n.boot}.}
#'    \item{n.neigh}{The number of waves carried out by the snowball sample.
#'          This is the same value from \code{sam.out$n.neigh}.}
#'    \item{seeds1}{A matrix of dimension \code{lenght(num.sam)} x \code{n.seeds} with
#'          the numeric seed ids. Each row corresponds to one LSMI. The rows are
#'          present in the same order as the ids in \code{num.sam}.}
#'    \item{nodes_of_LSMI}{A list of length \code{length(num.sam)} where each
#'          element is vector containing the numeric ids of the nodes sampled
#'          using the respective LSMI. The elements are present in the same
#'          order as the ids in \code{num.sam}.
#'          Note: nodes_of_LSMI is unreported when n.neigh equals zero.}
#' @export
#' @examples
#' net <- artificial_networks[[1]]
#' sam.out <- Oempdegreedistrib(net = net, n.seeds = 40, n.neigh = 1, num.sam = 1)
#' a <- bootdeg(sam.out = sam.out, n.boot = 50)

bootdeg <- function(sam.out, num.sam = sam.out$num.sam, n.boot = 1) {
      if (sam.out$n.neigh == 0) {
            # only information from the seeds
            res <- bootdeg0(sam.out, num.sam, n.boot)
      } else {
            res <- bootdegK(sam.out, num.sam, n.boot)
      }
      res
}
