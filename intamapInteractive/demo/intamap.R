data(intamap)
#observations = eurdepLoad(fname = "2006_1.csv",dataType = "eurdepAverage",sep = ",")
ck = unique(observations$isoCountry)


# For prediction locations
#data(boundaries)
#boundaries = boundaries[boundaries$COUNTRY != "RU",]
#projOrig = "+proj=laea +lat_0=48 +lon_0=9 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m"
#proj4string(boundaries) = CRS(projOrig)
#boundaries = spTransform(boundaries,CRS("+init=epsg:3035"))
#predictionLocations = spsample(boundaries,4000,"regular")
data(predictionLocations)
gridded(predictionLocations) = TRUE
data(countryBoundaries)

# set up intamap object:
krigingObject = list(
	pointData = observations,
	predictionLocations = predictionLocations,
  countryBoundaries = countryBoundaries,
	targetCRS = "+init=epsg:3035",
  formulaString = as.formula(value~1),
  ck = ck,
	params = getIntamapParams(isEmergency=TRUE)
)
class(krigingObject) = c("eurdep","automap")

# run test:
checkSetup(krigingObject)

# do interpolation steps:
# General preprocessing
krigingObject = preProcess(krigingObject)
krigingObject = estimateAnisotropy(krigingObject)
krigingObject = estimateParameters(krigingObject) # Does not yet take anisotropy into account
krigingObject = spatialPredict(krigingObject)  # Does not yet take anisotropy into account
krigingObject = postProcess(krigingObject)


# generate some output:
spplot(krigingObject$pointData,"value",col.regions = bpy.colors(),cuts=c(0,50,100,150,200,300,1000,10000,100000,1000000))
spplot(krigingObject$predictions,"var1.pred",col.regions = bpy.colors(),at=c(-2500,0,50,100,150,200,300,1000,10000,100000,1000000))
print(paste("Ratio :",krigingObject$anisPar$ratio," Direction",krigingObject$anisPar$direction))

