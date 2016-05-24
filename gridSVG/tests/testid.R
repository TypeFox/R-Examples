
library(grid)
library(gridSVG)

# A very simple test
dev.new(width=6, height=6)
# Some default settings
pushViewport(viewport(gp=gpar(col="black", fill=NA)))
grid.circle(r=4:1/8, name="circgrob")
grid.rect(name="rectgrob")
popViewport()

grid.export()
dev.off()
