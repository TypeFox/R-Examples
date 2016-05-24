#####################################
#                                   #
# RANK HISTOGRAM SIGNIFICANCE TESTS #
#                                   #
#####################################
TestRankhist <- function(rank.hist) {
#
# Conduct a series of significance tests for flatness of the rank histogram
#
# Usage: TestRankhist(rank.hist=rh)
#
# Arguments:
#   * rank.hist ... a vector of rank counts (see function `rankhist())`
#
# Return value:
#   * a dataframe with:
#     + rows    ... test.statistic, p-value
#     + columns ... pearson chi^2, jolliffe-primo slope, jolliffe-primo convex
#
# Author: 
#
#    Stefan Siegert 
#    s.siegert@exeter.ac.uk 
#    December 2013
#
# Example:
# 
#   ens <- matrix(rnorm(500), 100, 5)
#   obs <- rnorm(100)
#   rh <- rankhist(ens, obs)
#   rankhist.tests(rank.hist = rh)
#
# References: Pearson 1900 
#             Jolliffe & Primo 2008 http://dx.doi.org/10.1175/2007MWR2219.1
#

  o.i <- rank.hist
  N <- sum(o.i)
  J <- length(o.i) 
  i <- 1:J
  e.i <- N / J
  x.i <- (o.i - e.i) / sqrt(e.i)
  # pearson chi^2 test 
  X2 <- sum(x.i * x.i)
  p.chisq <- pchisq(X2, df=J-1, lower.tail=FALSE)
  # jolliffe-primo: 
  # linear contrast
  a <- 2 * sqrt(3 / (J^3 - J))
  b <- -(sqrt(3) * J + sqrt(3)) / sqrt(J * (J + 1) * (J - 1))
  x.lin <- a*i+b
  X2.lin <- sum(x.i * x.lin)^2 # should have chi^2(df=1)
  # squared contrast
  a <- 6 * sqrt(5 / (J^5 - 5 * J^3 + 4 * J))
  b <- -1 / 2 * (sqrt(5) * J^2 - sqrt(5)) / 
       (sqrt((J - 2) * (J - 1) * J * (J + 1) * (J + 2)))
  x.u <- a * (i - (J + 1) / 2)^2 + b
  X2.u <- sum(x.i * x.u)^2 # should have chisq(df=1)

  # return 
  ret.df <- data.frame(pearson.chi2=c(X2, p.chisq), jp.slope=c(X2.lin, pchisq(X2.lin, df=1, lower.tail=FALSE)), jp.convex=c(X2.u, pchisq(X2.u, df=1, lower.tail=FALSE)))
  rownames(ret.df) <- c("test.statistic", "p.value")
  ret.df
}


