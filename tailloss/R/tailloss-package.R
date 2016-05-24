#' Evaluate the Probability in the Upper Tail of the Aggregate Loss Distribution
#' 
#' Evaluate the probability in the upper tail of the aggregate loss distribution using different methods: Panjer recursion, Monte Carlo simulations, Markov bound, Cantelli bound, Moment bound, and Chernoff bound.

#' @details 
#' The package \code{tailloss} contains functions to estimate the exceedance probability curve of the aggregated losses. There are two `exact' approaches: Panjer recursion and Monte Carlo simulations, and four approaches producing upper bounds: the Markov bound, the Cantelli bound, the Moment bound, and the Chernoff bound. The upper bounds are useful and effective when the number of events in the catalogue is large, and there is interest in estimating the exceedance probabilities of exceptionally high losses. 

#' @author Isabella Gollini <isabella.gollini@@bristol.ac.uk>, and Jonathan Rougier.
#'
#' This work was supported by the Natural Environment Research Council [Consortium on Risk in the Environment: Diagnostics, Integration, Benchmarking, Learning and Elicitation (CREDIBLE); grant number NE/J017450/1]
#'
#' @references Gollini, I., and Rougier, J. C. (2015), "Rapidly bounding the exceedance probabilities of high aggregate losses", \url{http://arxiv.org/abs/1507.01853}.
#'
#' @examples
#' data(UShurricane)
#'
#' # Compress the table to millions of dollars
#'
#' USh.m <- compressELT(ELT(UShurricane), digits = -6)
#' s <- seq(1,40)

#' EPC <- matrix(NA, length(s), 6)
#' colnames(EPC) <- c("Panjer", "MonteCarlo", "Markov", 
#'  "Cantelli", "Moment", "Chernoff")
#' EPC[, 1] <- fPanjer(USh.m, s = s)[, 2]
#' EPC[, 2] <- fMonteCarlo(USh.m, s = s)[, 2]
#' EPC[, 3] <- fMarkov(USh.m, s = s)[, 2]
#' EPC[, 4] <- fCantelli(USh.m, s = s)[, 2]
#' EPC[, 5] <- fMoment(USh.m, s = s)[, 2]
#' EPC[, 6] <- fChernoff(USh.m, s = s)[, 2]

#' matplot(s, EPC, type = "l", lwd = 2, xlab = "s", ylim = c(0, 1), lty = 1:6,
#'   ylab = expression(plain(Pr)(S>=s)), main = "Exceedance Probability Curve")
#' zoombox(s, EPC, x0 = c(30, 40), y0 = c(0, .1), y1 = c(.3, .6), type = "l",
#'   lwd = 2, lty = 1:6)
#' legend("topright", legend = colnames(EPC), lty = 1:6, col = 1:6, lwd = 2)
#'

#' EPCcap <- matrix(NA, length(s), 6)
#' colnames(EPCcap) <- c("Panjer", "MonteCarlo", "Markov", 
#'  "Cantelli", "Moment", "Chernoff")
#' EPCcap[, 1] <- fPanjer(USh.m, s = s, theta = 2, cap = 5)[, 2]
#' EPCcap[, 2] <- fMonteCarlo(USh.m, s = s, theta = 2, cap = 5)[, 2]
#' EPCcap[, 3] <- fMarkov(USh.m, s = s, theta = 2, cap = 5)[, 2]
#' EPCcap[, 4] <- fCantelli(USh.m, s = s, theta = 2, cap = 5)[, 2]
#' EPCcap[, 5] <- fMoment(USh.m, s = s, theta = 2, cap = 5)[, 2]
#' EPCcap[, 6] <- fChernoff(USh.m, s = s, theta = 2, cap =  5)[, 2]

#' matplot(s, EPCcap, type = "l", lwd = 2, xlab = "s", ylim = c(0, 1), lty = 1:6,
#'   ylab = expression(plain(Pr)(S>=s)), main = "Exceedance Probability Curve")
#' zoombox(s, EPCcap, x0 = c(30, 40), y0 = c(0, .1), y1 = c(.3, .6), type = "l",
#'   lwd = 2, lty = 1:6)
#' legend("topright", legend = colnames(EPC), lty = 1:6, col = 1:6, lwd = 2)
#'

#'
#' @name tailloss-package
#' @aliases tailloss
#' @import MASS
#' @import graphics
#' @import stats
#' @docType package
NULL
