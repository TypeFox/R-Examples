# Calculate the AUC of ROC curve for ranked compounds
# x: a vector for scores
# y: a vector for labels
#
# e.g.)
# > x <- rnorm(100001) - 1:100001 * 0.00005
# > y <- c(rep(1,501), rep(0,length(x)-501))
# > auc(x, y, decreasing=TRUE)
#
# Ref.: Tom Fawcett, An introduction to ROC analysis. 
# Pattern Recognition Letters 27, 861-874 (2006)
#

auc <- function(x, y, decreasing=TRUE, top=1.0) {
  if ( length(x) != length(y) ){
    stop(paste("The number of scores must be equal to the number of labels."))
  }
  N <- length(y)
  n <- sum(y == 1)

  x_prev <- -Inf
  area <- 0
  fp = tp = fp_prev = tp_prev = 0
  ord <- order(x, decreasing=decreasing)
  for (i in seq_along(ord)) {
    j <- ord[i]
    if (x[j] != x_prev) {
      if( fp >= (N - n) * top ){
        rat <- ((N - n) * top - fp_prev) / (fp - fp_prev)
        area <- area + rat * (fp - fp_prev) * (tp + tp_prev) / 2
        return( area/( n *(N - n) * top) )
      }
      area <- area + (fp - fp_prev) * (tp + tp_prev) / 2
      x_prev <- x[j]
      fp_prev <- fp
      tp_prev <- tp
    }
    if (y[j] == 1) {
      tp <- tp + 1
    } else {
      fp <- fp + 1
    }
  }
  area <- area + (fp - fp_prev) * (tp + tp_prev) / 2
  return( area/(n*(N-n)) )
}
