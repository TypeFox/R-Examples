head(Wetsuits, 3)
dotPlot(~ Wetsuit, xlim = c(1.1, 1.8), cex = .25, data = Wetsuits) # to check for normality
dotPlot(~ NoWetsuit, xlim = c(1.1, 1.8), cex = .25, data = Wetsuits) # to check for normality

