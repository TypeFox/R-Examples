plotPR <-
function(piNet, ...){
  plot(c(0, piNet[, "Recall"]), c(0, piNet[, "Precision"]), type = "l",
       xlab = "recall", ylab = "precision", main = "PR Curve",
       xlim = 0:1, ylim = 0:1, ...)
}
