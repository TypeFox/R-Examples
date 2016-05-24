library(grid)
library(gridSVG)

x <- roundrectGrob(width=.8, height=.8, name="x", gp=gpar(lwd=5))

grid.newpage()
grid.draw(x)
grid.gradientFill("x", linearGradient(c("green", "yellow")))
grid.export("force-gradient.svg")

grid.newpage()
grid.draw(x)
grid.animate("x", "stroke-opacity"=1:0)
grid.export("force-animate.svg")

grid.newpage()
grid.draw(x)
grid.hyperlink("x", href="http://www.stat.auckland.ac.nz/")
grid.export("force-hyper.svg")

grid.newpage()
grid.draw(x)
grid.filter("x", filterEffect(feGaussianBlur(sd=5)))
grid.export("force-filter.svg")

grid.newpage()
grid.draw(x)
grid.patternFill("x", pattern(circleGrob()))
grid.export("force-pattern.svg")

grid.newpage()
grid.draw(x)
grid.garnish("x", title="tooltip")
grid.export("force-garnish.svg")
