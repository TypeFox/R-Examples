# Calculate a enrichment factor for the highly-ranked compounds
# x: a vector for scores
# y: a vector for labels
#
# e.g.)
# > x <- rnorm(100001) - 1:100001 * 0.00005
# > y <- c(rep(1,501), rep(0,length(x)-501))
# > enrichment_factor(x, y, top=0.05)

enrichment_factor <- function(x, y, top=0.05, decreasing=TRUE) {
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
    if(x[j] != x_prev) {
      if( fp + tp >= N * top ){
        n_right <- (fp - fp_prev) + (tp - tp_prev)
        rat <- (N * top - (fp_prev +  tp_prev) ) / n_right
        tp_r <- tp_prev + rat * (tp - tp_prev)
        return( (tp_r / (N*top) ) / (n/N) )
      }
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
#  n_right <- (fp - fp_prev) + (tp - tp_prev)
  return( 1 )
}

  
