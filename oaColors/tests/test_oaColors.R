library(oaColors)
library(grid)

plot(1:4, col = oaColors(color = c("red", "orange", "yellow", "green")), pch = 19, cex = 2)
lower = 0; upper = 8
grid.newpage()
plotvp <- viewport(x=0.1,
		y=0.1, xscale = c(lower, upper), yscale = c(lower, upper), 
		width = 0.85, height = 0.85, 
		name="plotRegion", 
		just = c("left", "bottom"))
pushViewport(plotvp)
grid.points(1:7, 1:7, pch = 19, gp = gpar(col = oaColors(c("pink", "blue", "yellow", "cyan", "green", "red", "orange" ), alpha = 0.8)   ) )
