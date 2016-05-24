library(grid)
library(gridSVG)

pdf(file = NULL)

pushViewport(viewport())
grid.rect()
grid.circle()
grid.text("hello, world!")
popViewport()

# All grobs and viewports should now be exported with a class attribute
# holding the value of their R class()
grid.export("test-classes.svg", addClasses = TRUE)
dev.off()
