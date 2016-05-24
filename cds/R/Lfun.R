#' Calculate Constrained Dual Scaling Loss
#' 
#' Calculate the loss function for constrained dual scaling.
#' 
#' @param a.cur The current value for a.
#' @param bkmat Current value of bkmat.
#' @param G Current value G.
#' @param Fr.cent Current value of the centred Fr.
#' @param n Number of respondents.
#' @param m Number of items.
#' @param q Number for rating scale categories so that the rating scale is \code{1:q}.
#' @param K Number of response style groups.
#' @param const Constant part of the loss function
#' @keywords multivariate
Lfun <- function(a.cur, bkmat, G, Fr.cent, n, m, q, const, K)	{
  mat <- matrix(a.cur, ncol = K, nrow = 2 * n)
	awts2 <- colSums(mat * mat * rbind(G, G))
	bwts2 <- colSums(bkmat * bkmat)
	last <- sum((mat * rbind(G, G))*(Fr.cent %*% bkmat))
  mq2 <- m + q - 2
	out <- const + 0.25 * mq2 * mq2 * sum(bwts2 * awts2) - mq2 * last
	return(out / const)
	}
