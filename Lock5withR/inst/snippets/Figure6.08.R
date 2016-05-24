plotFun(dnorm(x, 0, 1) ~x, x.lim = c(-4, 4), col = "black")
plotFun(dt(x, df = 15) ~x, add = TRUE, lty = 2)
plotFun(dt(x, df = 5) ~x, add = TRUE, lty = 3, col = "red")

