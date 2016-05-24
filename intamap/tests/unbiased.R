options(error = recover)
library(intamap)

data(meuse)
coordinates(meuse) = ~x+y
data(meuse.grid)
coordinates(meuse.grid) = ~x+y

meuse$zinc = log(meuse$zinc)

set.seed(112233)
krigingObject = createIntamapObject(
	observations = meuse,
	predictionLocations = spsample(meuse.grid,5,"regular"),
#	targetCRS = "+init=epsg:3035",
#	boundCRS = "+proj=laea +lat_0=48 +lon_0=9 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m",
#	boundCRS = boundCRS,
#	boundaries = boundaries,
  formulaString = as.formula("zinc~1"),
	params =  list(debug.level = 1),
  outputWhat = list(mean = TRUE, variance = TRUE, MOK=7,IWQSEL = 7,excprob = 7.0)
)
class(krigingObject) = c("automap")

checkSetup(krigingObject)
krigingObject = preProcess(krigingObject)
krigingObject = estimateParameters(krigingObject)
krigingObject = spatialPredict(krigingObject)
krigingObject = postProcess(krigingObject)
summary(krigingObject$outputTable)


class(krigingObject) = c("yamamoto")

checkSetup(krigingObject)
krigingObject = preProcess(krigingObject)
krigingObject = estimateParameters(krigingObject)
krigingObject = spatialPredict(krigingObject)
krigingObject = postProcess(krigingObject)
summary(krigingObject$outputTable)