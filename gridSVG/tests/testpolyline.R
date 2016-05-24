
library(grid)
require(gridSVG)

    grid.polyline(x=outer(c(0, .5, 1, .5), 5:1/5),
                  y=outer(c(.5, 1, .5, 0), 5:1/5),
                  id.lengths=rep(4, 5),
                  gp=gpar(col=1:5, lwd=3))
  
grid.export("polyline.svg")
