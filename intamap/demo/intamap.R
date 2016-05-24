data(intamap)
require(rworldmap)
countryBoundaries = getMap()
proj4string(observations) = "+init=epsg:4236"
observations$country = countryBoundaries$ISO3[overlay(observations, countryBoundaries)]
# removing observations that have fallen outside country borders, mostly
# because of the use of a low resolution map
observations = observations[!is.na(observations$country),]
# Keeping only boundaries of countries with observations
countryBoundaries = countryBoundaries[countryBoundaries$ISO3 %in% observations$country,]


predictionLocations = spsample(countryBoundaries,4000,"regular")

# set up intamap object:
krigingObject = createIntamapObject(
	observations = observations,
	predictionLocations = predictionLocations,
  countryBoundaries = countryBoundaries,
	targetCRS = "+init=epsg:3035",
  intCRS = "+init=epsg:3035",
  formulaString = as.formula(obs~1),
  class = "automap",
  params = list(confProj = TRUE)
)

# run test:
checkSetup(krigingObject)

# do interpolation steps:
# General preprocessing
krigingObject = preProcess(krigingObject)
krigingObject = estimateAnisotropy(krigingObject)
krigingObject = estimateParameters(krigingObject) 
krigingObject = spatialPredict(krigingObject)  
krigingObject = postProcess(krigingObject)


# generate some output:
spplot(krigingObject$observations,"obs",col.regions = bpy.colors())
spplot(krigingObject$predictions,"var1.pred",col.regions = bpy.colors())
print(paste("Ratio :",krigingObject$anisPar$ratio," Direction",krigingObject$anisPar$direction))

