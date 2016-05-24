#' Calculate Loss for G Update
#' 
#' Loss function used for updating G. This is not equivalent to the original loss function, as 
#' only a part of the total loss depends on G.
#' 
#' @param G The current value for G.
#' @param a.cur The current value for a.
#' @param bwts2 The current value of the squared b weights.
#' @param Fr.bk Current product between Fr.cent and bk. 
#' @param n Number of respondents.
#' @param m Number of items.
#' @param q Number for rating scale categories so that the rating scale is \code{1:q}.
#' @param K Number of response style groups.
#' @keywords multivariate
Lfun.G.upd <- function(G, a.cur, bwts2, Fr.bk, n, m, q, K) { #, i = NULL) #awts2.base = NULL
  mat <- matrix(a.cur, nrow = 2 * n, ncol = K)
  awts2 <- colSums(mat * mat * rbind(G, G))
	last <- colSums(mat * rbind(G, G)* Fr.bk)
	outvec <- 0.25*(m + q - 2)* bwts2 * awts2 - last
  out <- sum(outvec)
#   rm(awts2, last, mat)
	list(out = out, kloss = outvec)
}
