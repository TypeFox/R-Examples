plot.gSignature <-
function(x, ...) {
  tmp <- survfit(stData ~ x$classification)
  mmain <- paste("signature:", x$signatureName)
  xxlab <- paste("tValue(Log-Rank test) = ", round(x$tValue, 2), ", log(pValue, 10) = ", round(log(x$pValue, 10), 3), sep = "")
  plot(tmp, main = mmain, xlab = xxlab, col = c("green", "red"))
#  NextMethod("plot")
}
