plotFun(dt(x, df = 15) ~x, x.lim = c(-4, 4))
plotDist("t", params = list(df=15), type = c("h","l"), groups =  x > 1.5, lty = 1)
ladd(grid.text("1.5",1.5,.2, default.units = "native", hjust = 0))

