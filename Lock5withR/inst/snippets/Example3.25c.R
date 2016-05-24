cdata( ~ diffmean, 0.95, data = BootE)
dotPlot(~ diffmean, width = .25, cex = .75, groups =(-1.717 <= M & M <= 7.633), 
        xlab = "Difference in mean", data = BootE)

