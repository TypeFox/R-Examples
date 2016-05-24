P.05 <- do(2000) * rflip(50, .05)
dotPlot(~prop, width = .02, cex = 25, data = P.05)

P.10 <- do(2000) * rflip(50, .10)
dotPlot(~prop, width = .02, cex = 15, data = P.10)

P.25 <- do(2000) * rflip(50, .25)
dotPlot(~prop, width = .02, cex = 10, data = P.25)

P.50 <- do(2000) * rflip(50, .50)
dotPlot(~prop, width = .02, cex = 5, data = P.50)

P.90 <- do(2000) * rflip(50, .90)
dotPlot(~prop, width = .02, cex = 10, data = P.90)

P.99 <- do(2000) * rflip(50, .99)
dotPlot(~prop, width = .02, cex = 25, data = P.99)

