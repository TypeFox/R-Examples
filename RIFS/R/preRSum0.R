# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function:
# preRSum0() - a sample of fractal points in R^n 
#              generated with a matrix of random sums.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# n - a sample size;
# mu - a partition coefficient for iterative segments;
# eps - an error of a random sum of a numerical series;
# Z - a matrix of protofractal points coordinates;
# P - a probability distribution of protofractal points.
# Variables:
# k - a number of protofractal points; S - a matrix of random sums;
# m, X - a length and elements of the numerical series; 
# l - a sample of protofractal indexes.
# Value:
# list(pre, proto, distr)
# $pre - a sample of fractal points;
# $proto - a set of protofractal points;
# $distr - a distribution of probabilities & coefficients.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
preRSum0 <- function(n=10000, mu=1, eps=1e-9, Z=R2ngon(), 
                     P=rep(1/nrow(Z), times=nrow(Z))) {
  k <- nrow(Z)
  S <- array(0, dim=c(n,k))
  m <- 1 - log(eps*mu)/log(1+mu)
  X <- mu/(mu+1)^seq(m)
  for (i in seq(n)) {
    l <- sample.int(k, size=m, prob=P, replace=TRUE)
    for (j in seq(k)) 
      S[i,j] <- sum(X[l==j])
  }
  return(list(pre=S%*%Z, proto=Z,
              distr=cbind(p=P, mu=rep(mu, times=k))))
}