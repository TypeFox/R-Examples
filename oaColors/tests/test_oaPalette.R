library(oaColors)
library(grid)

n = 7
myColors <- oaPalette(n, alpha = 0.7)


lower = 0; upper = n+1
grid.newpage()
plotvp <- viewport(x=0.1,
		y=0.1, xscale = c(lower, upper), yscale = c(lower, upper), 
		width = 0.85, height = 0.85, 
		name="plotRegion", 
		just = c("left", "bottom"))
pushViewport(plotvp)
grid.points(1:n, 1:n, pch = 19, gp = gpar(col = myColors  ) )