plot.outlierProbs <- function(x,...) {
  probs <- x$outlier.prob
  names(probs) <- x$slab
   dotchart(rev(probs), xlim = c(0,1),xlab="Outlier Probability")
  abline(v=0.0,lty=3)
  abline(v=1.0,lty=3)
}