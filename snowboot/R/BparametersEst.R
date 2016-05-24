#' Summary of the Bootstrap Degree Distribution
#'
#' This function provides summary statistics of a bootstrap degree
#' distribution.
#' @param outBootdeg A list that is the output of \code{\link{bootdeg}}
#' @references Efron, B. (1979). Bootstrap methods: another look at the
#'  jackknife. The annals of Statistics, 1-26.
#' @references Thompson, M. E., Ramirez Ramirez, L. L., Lyubchich, V. and
#'  Gel, Y. R. (2015), Using the bootstrap for statistical inference
#'  on random graphs. Can J Statistics. doi: 10.1002/cjs.11271
#' @return A list consisting of:
#'    \item{mean}{An array of dimension
#'    \code{c(length(outBootdeg$num.sam),outBootdeg$n.boot,3)}.
#'    The last dimension, of 3, is for the three different methods of obtaining
#'    the empirical degree distribution from \code{outBootdeg$empd}
#'    (see output empd from \code{\link{bootdeg}} for details).
#'    The (i,j,k)-th element in the array is an estimate of mean degree for the
#'    i-th LSMI sample, j-th bootstrap replication, and k-th empirical distribution
#'    from \code{outBootdeg$empd}.}
#'    \item{quartiles}{An array of dimension
#'    \code{c(length(outBootdeg$num.sam), 3, outBootdeg$n.boot, 3)}.
#'    The last dimension, of 3, is for the three different methods of estimation
#'    from \code{outBootdeg$empd} (see output empd from \code{\link{bootdeg}}
#'    for details). The second dimension, of 3, corresponds to
#'    the quartiles (.25, .5, .75). The (i,j,k,l)-th element in the array is
#'    an estimate of j-th quartile for the i-th LSMI sample,
#'    k-th bootstrap replication, and l-th empirical distribution from
#'    \code{outBootdeg$empd}.}
#'    \item{rfreq}{An array of dimension
#'    \code{c(length(outBootdeg$num.sam), 5, outBootdeg$n.boot, 3)}.
#'    The last dimension, of 3, is for the three different methods of estimation
#'    from \code{outBootdeg$empd}
#'    (see output empd from \code{\link{bootdeg}} for details.).
#'    The second dimension, of 5, corresponds to degree values: 0, 1, 2, 3, 4.
#'    The (i,j,k,l)-th element in the array is the proportion of nodes
#'    with degree j in the i-th LSMI sample, k-th bootstrap replication,
#'    and l-th empirical distribution from \code{outBootdeg$empd}.}
#'    \item{deciles}{An array of dimension \code{c(length(outBootdeg$num.sam), 9,
#'          outBootdeg$n.boot, 3)}. The last dimension, of 3, is for the three
#'          different methods of estimation from \code{outBootdeg$empd}
#'          (see output empd from \code{\link{bootdeg}} for details.). The second
#'          dimension, of 9, corresponds to the deciles (.1, .2, ... , .9).
#'          The (i,j,k,l)-th element in the array is an estimate of j-th
#'          decile for the i-th LSMI sample, k-th bootstrap replication,
#'          and l-th empirical distribution from \code{outBootdeg$empd}.}
#'    \item{num.sam}{Numeric indices corresponding to LSMI samples used for bootstrap.
#'          See value \code{num.sam} from \code{\link{bootdeg}}.}
#'    \item{seeds1}{A matrix of dimension \code{length(num.sam)} x \code{n.seeds} with
#'          the numeric seed ids. Each row corresponds to one LSMI. The rows are
#'          present in the same order as the ids in \code{num.sam}.
#'          See value \code{seeds1} from \code{\link{bootdeg}}.}
#'    \item{nodes_of_LSMI}{A list of length \code{length(num.sam)} where each
#'          element is vector containing the numeric ids of the nodes sampled
#'          using the respective LSMI. The elements are present in the same
#'          order as the ids in \code{num.sam}.
#'          Note: nodes_of_LSMI is unreported when n.neigh equals zero.
#'          See value \code{nodes_of_LSMI} from \code{\link{bootdeg}}.}
#' @export
#' @examples
#' net <- artificial_networks[[1]]
#' sam.out <- Oempdegreedistrib(net = net, n.seeds = 40, n.neigh = 1, num.sam = 1)
#' outBootdeg <- bootdeg(sam.out = sam.out, n.boot = 50)
#' a <- BparametersEst(outBootdeg)

BparametersEst <- function(outBootdeg) {
      # outBootdeg is output of bootdeg
      tn.sam <- length(outBootdeg$num.sam)
      if (outBootdeg$n.neigh == 0) {
            n.dist <- 1  #n.dist is the number of different emp distr.
      } else {
            n.dist <- 3
      }

      mean <- array(NA, c(tn.sam, outBootdeg$n.boot, n.dist))
      quartiles <- array(NA, c(tn.sam, 3, outBootdeg$n.boot, n.dist))
      rfreq <- array(NA, c(tn.sam, 5, outBootdeg$n.boot, n.dist))  #freq of 0,1,2,3,4
      deciles <- array(NA, c(tn.sam, 9, outBootdeg$n.boot, n.dist))  #10,20,30,40,50,60,70,80,90

      for (m in 1:tn.sam) {
            w <- 1
            in.while <- TRUE
            while (in.while) {
                  if (outBootdeg$n.neigh == 0) {
                        empd <- outBootdeg$empd[[m]]$empd.seeds
                        in.while <- FALSE
                  } else {
                        if (w == 1) {
                              empd <- outBootdeg$empd[[m]]$empd.w.p0s
                        } else if (w == 2) {
                              empd <- outBootdeg$empd[[m]]$empd.nw.p0sEkb
                        } else if (w == 3) {
                              empd <- outBootdeg$empd[[m]]$empd.nw.p0sEks
                              in.while <- FALSE
                        }
                  }
                  vals <- outBootdeg$values[[m]]  #as.numeric(colnames(empd))  ###checar

                  if (dim(empd)[2] != length(as.vector(vals))) {
                        # browser()
                        mean[m, , w] <- t(empd) %*% as.vector(vals)
                  } else mean[m, , w] <- empd %*% as.vector(vals)

                  cempd <- t(apply(empd, 1, FUN = cumsum))
                  quartiles[m, , , w] <- t(sapply(X = c(0.25, 0.5, 0.75), FUN = cempdpercentile, cempd = cempd, vals = vals))
                  rfreq[m, , , w] <- t(sapply(X = 0:4, FUN = distribvalsmat, empd = empd, vals = vals))
                  deciles[m, , , w] <- t(sapply(X = seq(0.1, 0.9, by = 0.1), FUN = cempdpercentile, cempd = cempd, vals = vals))
                  w <- w + 1
            }  #while(w<=3)
      }  #for (m in 1:...)
      # browser()
      list(mean = mean, quartiles = quartiles, rfreq = rfreq, deciles = deciles,
           n.dist = n.dist, num.sam = outBootdeg$num.sam,
           seeds1 = outBootdeg$seeds1, nodes_of_LSMI = outBootdeg$nodes_of_LSMI)
}  #function
