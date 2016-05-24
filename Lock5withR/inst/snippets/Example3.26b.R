cdata( ~ cor, 0.98, data = BootM)
dotPlot( ~ cor, width = .005, groups = (-.940 <= cor & cor <= -.705), 
         xlab = "r", data = BootM)

