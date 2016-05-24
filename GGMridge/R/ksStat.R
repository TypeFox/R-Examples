#----------------------------------------------------------------------#
# Calculates Kolmogorov-Smirnov statistic for p-values                 #
# Compare p values with uniform distribution                           #
# Modified from ks.test function                                       #
#----------------------------------------------------------------------#
# Input :                                                              #
#                                                                      #
#  p  A numeric vector indicating pvalues.                             #
#                                                                      #
# Output :                                                             #
#                                                                      #
#  Kolmogorov-Smirnov statistic                                        #
#----------------------------------------------------------------------#
ksStat <- function(p) {

  p <- p[!is.na(p)]
  n <- length(p)

  p <- punif(q = sort(p)) - {0L:{n - 1L}} / n

  ksStat <-  max( c(p, 1.0/n - p))

  return(ksStat)
}
