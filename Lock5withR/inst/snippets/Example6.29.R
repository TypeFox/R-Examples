head(Wetsuits, 3)
t.test(Wetsuits$Wetsuit, Wetsuits$NoWetsuit, paired = TRUE)
dotPlot(Wetsuits$Wetsuit - Wetsuits$NoWetsuit, width = .01, cex = .3)

