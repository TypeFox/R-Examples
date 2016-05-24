CGaPlotDiag <-
function(M, variable = "lambda", pos = 1) {
  if (variable != "lambda" && variable != "u" && variable != "c"
      && variable != "epsilon" && variable != "theta") {
    stop ("Error: 'variable' must be either 'lambda', 'u', 'c',
          'epsilon' or 'theta'.")
  }
  tol = .Machine$double.eps ^ 0.5
  MAT <- M$summary
  p <- M$p
  K <- M$K
  a <- 0
  if (length(MAT[, 1]) == 3 * K - 1 + p) {
    a <- 1
  }
  if (pos < 0 || pos > K || abs(pos - round(pos)) > tol ) {
    stop ("Invalid position.")
  }
  if (pos > (K - 1) && (variable == "u" || variable == "c")) {
    stop ("Invalid position.")
  }
  if (pos > p && variable == "theta"){
    stop ("Invalid position.")
  }
  if (variable == "epsilon" && pos != 1) {
    warning("'epsilon' has only one entry (1). Graphics shown for epsilon_1.")
    pos <- 1
  }
  if (variable == "lambda") {
    b <- 0
  }
  if (variable == "u") {
    b <- K
    if (mean(MAT[pos + b, ]) == 0) {
      stop ("Plots for 'u' are not available.")
    }
  }
  if (variable == "c") {
    b <- 2 * K - 1
    if (mean(MAT[pos + b, ]) == 0) {
      stop ("Plots for 'c' are not available.")
    }
  }
  if (variable == "epsilon") {
    b <- 3 * K - 2
  }
  if (variable == "theta") {
    b <- 3 * K - 2 + a
  }
  if (var(MAT[pos + b, ]) == 0) {
    stop ("Plots are not available.")
  }  
  par(mfrow = c(2, 2))
  ## Trace
  plot(MAT[pos + b, ], type = "l", xlab = "Iteration", ylab = "", 
       col = "slateblue4")
  mtext("Trace", line = 1, ps = 2, cex = 1, font = 1)
  ## Ergodic Mean
  p.erg <- cumsum(MAT[pos + b, ]) / 1:length(MAT[1, ])
  plot(p.erg, type = "l", xlab = "Iteration", ylab = "", col = "slateblue4")
  mtext("Ergodic mean", line = 1, ps = 2, cex = 1, font = 1)
  ## Autocorrelation function
  acf(MAT[pos + b, ], main = "", ylab = "", lwd = 2)
  mtext("Autocorrelation function", line = 1, ps = 2, cex = 1, font = 1)
  ## Histogram
  hist(MAT[pos + b, ], main = "", xlab = "", col = "lightblue", freq = FALSE)
  mtext("Histogram", line = 1, ps = 2, cex = 1, font = 1)
  par(mfrow = c(1, 1))
  mtext(paste(variable, "_", pos), line = 2.5, ps = 2, cex = 1.25, font = 2)
}
