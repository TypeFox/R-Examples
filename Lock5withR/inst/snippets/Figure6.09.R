plotFun(dt(x, df = 15) ~x, x.lim = c(-4, 4))
plotDist("t", params = list(df=15), type = c("h","l"), groups = (-2.131 < x & x < 2.131), lty = 1)
ladd(grid.text("2.131",2.1,.1, default.units = "native", hjust = 0))

