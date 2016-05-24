library(grid)
library(gridSVG)

dev.new(width=6, height=6)
grid.rect(name = "mainrect")
grid.comment("This is a comment", "mainrect")
grid.export("comment-test.svg")
dev.off()
