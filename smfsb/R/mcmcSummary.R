
mcmcSummary <- function (mat, rows = 4, lag.max=100, bins=30, show = TRUE, plot = TRUE) 
{
  d = dim(mat)
  p = d[2]
  summ=summary(mat)
  if (show==TRUE) {
    message(paste("N =", d[1], "iterations"))
    print(summ)
    message("Standard deviations:")
    print(apply(mat, 2, sd))
  }
  if (plot==TRUE) {
    names = colnames(mat)
    op = par(mfrow = c(rows, 3))
    for (i in 1:p) {
      plot(ts(mat[, i]), main = names[i], ylab = "Value", xlab = "Iteration")
      acf(mat[, i], lag.max = lag.max, main = names[i], ci=0, ylim=c(0,1))
      hist(mat[, i], bins, main = names[i], xlab = "Value", freq = FALSE)
    }
    par(op)
  }
  invisible(summ)
}

# eof
