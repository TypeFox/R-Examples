library(grid)
library(gridSVG)

pdf(file = NULL)

pushViewport(viewport())
grid.rect()
grid.circle(name = "mycircle")
grid.text("hello, world!")
popViewport()

# All IDs should now be exported with the prefix "TESTPREFIX"
grid.export("test-prefixes.svg", prefix = "TESTPREFIX")
dev.off()
