plotFun(dt(x, df = 49) ~x, x.lim = c(-4, 4))
plotDist("t", params = list(df=49), type = c("h","l"), groups = (-3.14 < x & x < 3.14), lty = 1)
ladd(grid.text("3.14",3,.05, default.units = "native", hjust = 0))

