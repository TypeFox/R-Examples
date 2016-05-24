n10 <- do(2000) * rflip(10, .10)
dotPlot(~prop, width = .1, cex = 25, data = n10)

n25 <- do(2000) * rflip(25, .10)
dotPlot(~prop, width = .04, cex = 10, data = n25)

n200 <- do(2000) * rflip(200, .10)
dotPlot(~prop, width = .005, cex = 5, data = n200)

