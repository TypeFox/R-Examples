## ----loading, results='hide', echo=FALSE, warning=FALSE, message=FALSE----
library(SOMbrero)

## ------------------------------------------------------------------------
my.first.grid <- initGrid(dimension=c(5,6), topo="square", dist.type="maximum")

## ----plotFirstGrid-------------------------------------------------------
print(my.first.grid)

## ------------------------------------------------------------------------
summary(my.first.grid)

## ----plotRainbowGrid-----------------------------------------------------
# plot without any color specification: squares are filled in white color
plot(my.first.grid)

# generating colors from the rainbow() function
my.colors <- rainbow(5*6)
plot(my.first.grid, neuron.col=my.colors)

