#' Lower bound estimator for discrete powerlaw distributions
#'
#' Lower bound estimator for discrete powerlaw distributions based on the distances between probability mass functions.
#' @param o Discrete powerlaw object.
#' @param g A guess on the true value of the lower bound.
#' @param c Confidence on the guess. A value between 1 and 100. By default is set to 90.
#' @param k Number of computations after a local minimum in the KS statistics is reached.
#' @param xmax Max value considered in the estimation of the lower bound.
#' @keywords discrete powerlaw lower bound estimator
#' @references A. Bessi, Speeding up lower bound estimation in powerlaw distributions, arXiv
#' @export getXmin2
#' @examples
#' x = moby
#' o = displo(x)
#' est = getXmin2(o)

getXmin2 = function(o, g = 1, c = 90, k = 5, xmax = 1e5)
{
  est = list()

  x = o$x
  N = o$nx
  g = g - (g * (100-c) / 100)

  # initialize
  xmins = o$ux[o$ux<=xmax]
  START = xmins[which.min(abs(xmins - g))]
  if (START > g){
    START = xmins[which.min(abs(xmins - g)) - 1]
  }
  if (length(START) == 0)
  {
    START = 1
  }
  xmins = xmins[which(xmins==START):length(xmins)]
  L = length(xmins)
  KS = numeric()
  alpha = numeric()
  xu = sort(x)
  len_xu = length(xu)

  for (i in 1:L)
  {
    n = length(xu[xu>=xmins[i]])
    q = xu[(N-n+1):len_xu]
    q = q[q <= xmax]

    # KS
    S = pmf(q)$y
    alpha = c(alpha, 1 + length(q) / sum(log(q/(xmins[i]-0.5))))
    P = ddispl(unique(q), xmin = xmins[i], alpha = alpha[i])
    KS = c(KS, max(abs(P-S)))
    if (i >= (k + 1))
      if( length(which((diff(KS[(i-k):i])>0) == TRUE)) == k )
        break

  }
  est$xmin = xmins[which.min(KS)]
  est$alpha = alpha[which.min(KS)]

  o$xmin = xmins[which.min(KS)]
  o$alpha = alpha[which.min(KS)]
  o$sigma = (alpha[which.min(KS)] - 1) / sqrt(N)

  return(est)
}
