library(hexbin)

## a variation on Nicholas' post to bioconductor & example(hexViewport)
set.seed(545)
x <- rnorm(2^15)
y <- 3*x - .2*x^2 + rnorm(2^15)
hbin <- hexbin(x,y)

##
hp <- hexViewport(hbin, newpage = TRUE)
pushHexport(hp)
library("grid")
grid.rect()
grid.xaxis()
grid.yaxis()
grid.hexagons(hbin, style = "centroid")
hloess <- loess(y ~ x, data = hcell2xy(hbin), weights = hbin @ count)
xx <- seq(hbin@xbnds[1], hbin@xbnds[2], length = 500)
grid.lines(xx, predict(hloess, xx),
           gp = gpar(col = 'red', lwd = 2), default.units = "native")
popViewport()
