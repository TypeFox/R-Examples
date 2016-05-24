set.seed(15331)
library(intamap)
data(meuse)
data(meuse.grid)
coordinates(meuse) = ~x+y
coordinates(meuse.grid) = ~x+y

meuse$value=meuse$zinc
output = interpolate(meuse, meuse.grid, list(mean=T, variance=T),methodName = "transGaussian")
summary(t(output$outputTable))


