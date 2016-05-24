plotROC <-
function(piNet, ...){
  plot(c(0, piNet[, "FPR"], 1), c(0, piNet[, "Recall"], 1), type = "l",
       xlab = "FP rate", ylab = "TP rate", main = "ROC Curve",
       xlim = 0:1, ylim = 0:1, ...)
  lines(0:1, 0:1, col = "black")
}
