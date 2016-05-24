library(intamap)

# observations = Something from Java...
# Until then we use the Meuse data:
data(meuse)
observations = data.frame(x = meuse$x,y = meuse$y,value = log(meuse$zinc))
# If you send a field just with 3 columns (x,y & z), we can let R figure
# out itself which names they have, for creation of a spatial object:
obsNames = names(observations)
coordinates(observations) = as.formula(paste("~",obsNames[1], "+", obsNames[2]))
set.seed(13531)
predictionLocations = spsample(observations, 50, "regular")
proj4string(observations) <- CRS("+proj=stere +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.999908 +x_0=155000 +y_0=463000 +ellps=bessel +units=m")
proj4string(predictionLocations) <- CRS("+proj=stere +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.999908 +x_0=155000 +y_0=463000 +ellps=bessel +units=m")

# We dont know the projection of the data at this stage, assume it is
# somehow metric

krigingObject = createIntamapObject(
	observations = observations,
	predictionLocations = predictionLocations,
  targetCRS = "+init=epsg:3035",
#	boundCRS = "+proj=laea +lat_0=48 +lon_0=9 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m",
#	boundCRS = boundCRS,
#	boundaries = boundaries,
  formulaString = as.formula(paste(obsNames[3],"~1")),
	params =  list(confProj = TRUE, thresh = quantile(observations$value,0.9)),
  outputWhat = list(mean=TRUE, variance=TRUE, excprob = 5.9, cumdistr = 5.9, 
		quantile = .1), class = "automap"

)

checkSetup(krigingObject)
krigingObject = preProcess(krigingObject)
krigingObject = estimateParameters(krigingObject)
krigingObject = spatialPredict(krigingObject)
krigingObject = postProcess(krigingObject)

summary(krigingObject$outputTable)

