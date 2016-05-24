library(intamap)

set.seed(13131)

# set up data:
data(meuse)
coordinates(meuse) = ~x+y
meuse$value = log(meuse$zinc)
data(meuse.grid)
gridded(meuse.grid) = ~x+y
proj4string(meuse) = CRS("+init=epsg:28992")
proj4string(meuse.grid) = CRS("+init=epsg:28992")

# set up intamap object:
idwObject = createIntamapObject(
	observations = meuse,
	formulaString=as.formula(zinc~1),
  predictionLocations = meuse.grid,
	targetCRS = "+init=epsg:3035",
	constantBias = 0,
	classes = "idw"
)

# run test:
checkSetup(idwObject)

# do interpolation steps:
idwObject = preProcess(idwObject)
idwObject = estimateParameters(idwObject, idpRange = seq(0.25,2.75,.5), nfold=3) # faster
idwObject = spatialPredict(idwObject)
idwObject = postProcess(idwObject)

# generate some output:
summary(as.data.frame(idwObject$outputTable))
