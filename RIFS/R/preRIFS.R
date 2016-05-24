# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Function:
# preRIFS() - a sample of fractal points in R^n 
#             generated with a random iterated function system.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Arguments:
# n - a sample size;
# Z - a matrix of protofractal points coordinates;
# P - a probability distribution of protofractal points;
# M - a partition coefficients distribution of protofractal points.
# Variables:
# z - a random sample of indexes for protofractal points; 
# X - a random sample of fractal points (prefractal points).
# Value:
# list(pre, proto, distr, index)
# $pre - a sample of fractal points;
# $proto - a set of protofractal points;
# $distr - a distribution of probabilities & coefficients;
# $index - a sample of protofractal indexes.
# # # # # # # # # # # # # # # # # # # # #
preRIFS <- function(n=10000, Z=R2ngon(), 
                    P=rep(1/nrow(Z), times=nrow(Z)),
                    M=rep(1, times=nrow(Z))) {
  z <- sample.int(nrow(Z), size=n, prob=P, replace=TRUE)
  X <- array(0, dim=c(n, ncol(Z)))
  X[1,] <- colMeans(Z)
  for (i in 2: n)
    X[i,] <- (X[i-1,] + M[z[i]]*Z[z[i],])/(1 + M[z[i]])
  return(list(pre=X, proto=Z,
              distr=cbind(p=P, mu=M),
              index=z))
}