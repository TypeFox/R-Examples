plot.ICA <-
function(x, ..., xlab = "Iteration", ylab = "Value",
                     main = "ICA History", col = "red") {
  plot(x$trace, xlab = xlab, ylab = ylab, main = main, type = "l", col = col)
}
