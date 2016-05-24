library(grid)
library(gridSVG)
dev.new(width=6, height=6)
# Some default settings
pushViewport(viewport(gp=gpar(col="black", fill=NA)))

# A plot with overlapping polygons where the overlap needs to be
# coloured differently

# Implemented using transparency

# Some dummy data
set.seed(1000)
x <- sort(runif(20, 1, 10))
y1 <- (x - 5) + rnorm(20)
y2 <- -(x - 5) + rnorm(20)

# Some "analysis"
lm1 <- lm(y1 ~ x)
lm2 <- lm(y2 ~ x)

# Some calculated values to plot
p1 <- predict(lm1, interval="confidence", type="response")
p2 <- predict(lm2, interval="confidence", type="response")

pushViewport(plotViewport(c(5, 5, 4, 2)))
pushViewport(dataViewport(x, c(p1[,2], p2[,2], p1[,3], p2[,3])))
grid.rect()
grid.xaxis()
grid.yaxis()
grid.points(x, y1)
grid.points(x, y2, pch=3)
grid.lines(x, p1[,1], default.units="native")
grid.lines(x, p1[,2], default.units="native")
grid.lines(x, p1[,3], default.units="native")
grid.lines(x, p2[,1], default.units="native")
grid.lines(x, p2[,2], default.units="native")
grid.lines(x, p2[,3], default.units="native")
# overlapping polygons
grid.polygon(c(x, rev(x)), c(p1[,2], rev(p1[,3])),
             gp=gpar(fill="red", alpha=0.5),
             default.units="native")
grid.polygon(c(x, rev(x)), c(p2[,2], rev(p2[,3])),
             gp=gpar(fill="green", alpha=0.5),
             default.units="native")


popViewport(2)

popViewport()

grid.export("transparency.svg")
dev.off()
