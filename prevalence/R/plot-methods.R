plot.betaPERT <-
function(x, y, ...){
  main <-
    paste("Beta4(", round(x$alpha, 3), ", ", round(x$beta, 3), ", ",
          round(x$a, 3), ", ", round(x$b, 3), ")", sep = "")
  x_val <-
    seq(x$a, x$b, length.out = 1000)
  y_val <-
    with(x, dbeta(seq(0, 1, length.out = 1000), alpha, beta) * (b - a) + a)

  plot(x_val, y_val, col = "blue", type = "l", lwd = 2,
       xlab = "x", ylab = "density", main = main)
}

plot.betaExpert <-
function(x, y, ...){
  main <-
    paste("Beta(", round(x$alpha, 3), ", ", round(x$beta, 3), ")", sep = "")
  with(x,
       curve(dbeta(x, alpha, beta),
             col = "blue", lwd = "2",
             xlab = "x", ylab = "density", main = main))
}