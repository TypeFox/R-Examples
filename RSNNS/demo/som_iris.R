library(RSNNS)

#data(snnsData)
#inputs <- snnsData$som_cube.pat

data(iris)
inputs <- normalizeData(iris[,1:4], "norm")

model <- som(inputs, mapX=16, mapY=16, maxit=500,  calculateActMaps=TRUE, targets=iris[,5])

par(mfrow=c(3,3))

for(i in 1:ncol(inputs)) plotActMap(model$componentMaps[[i]], col=rev(topo.colors(12)))

plotActMap(model$map, col=rev(heat.colors(12)))
plotActMap(log(model$map+1), col=rev(heat.colors(12)))
persp(1:model$archParams$mapX, 1:model$archParams$mapY, log(model$map+1), 
                          theta = 30, phi = 30, expand = 0.5, col = "lightblue")

plotActMap(model$labeledMap)

model$componentMaps
model$labeledUnits
model$map

names(model)
