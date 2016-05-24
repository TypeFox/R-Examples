# CompareHypergeo.R
# This demo shows the difference between the three distributions:
# 1. Wallenius' noncentral hypergeometric distribution
# 2. Fisher's noncentral hypergeometric distribution
# 3. The (central) hypergeometric distribution

require(BiasedUrn)
require(stats)

ComparePlot <- function(m1, m2, n, odds) {
  xmin <- minHypergeo(m1, m2, n)
  xmax <- maxHypergeo(m1, m2, n)
  x <- xmin : xmax
  wnc <- dWNCHypergeo(x, m1, m2, n, odds)
  fnc <- dFNCHypergeo(x, m1, m2, n, odds)
  hyp <- dhyper(x, m1, m2, n)
  plot   (x, wnc, type="l", col="blue", 
     main = "Hypergeometric distributions", 
     sub = "Blue = Wallenius, Red = Fisher, Green = Central",
     xlab = "x", ylab = "Probability")
  points (x, fnc, type="l", col="red")
  points (x, hyp, type="l", col="green")
}

ComparePlot(80, 60, 100, 0.5)
