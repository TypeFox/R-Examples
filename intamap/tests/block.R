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
predictionLocations = spsample(observations, 10, "regular")
gridded(predictionLocations) = TRUE
cellsize = predictionLocations@grid@cellsize
cs = predictionLocations@grid@cellsize[1]/2

# We dont know the projection of the data at this stage, assume it is
# somehow metric

Srl = list()
for (i in 1:dim(coordinates(predictionLocations))[1]) {
  pt1 = coordinates(predictionLocations)[i,]
  x1 = pt1[1]-cs
  x2 = pt1[1]+cs
  y1 = pt1[2]-cs
  y2 = pt1[2]+cs

  boun = data.frame(x=c(x1,x2,x2,x1,x1),y=c(y1,y1,y2,y2,y1))
  coordinates(boun) = ~x+y
  boun = Polygon(boun)
  Srl[[i]] = Polygons(list(boun),ID = as.character(i))
}
predictionLocations = SpatialPolygons(Srl)


krigingObject = createIntamapObject(
	observations = observations,
	predictionLocations = predictionLocations,
#	targetCRS = "+init=epsg:3035",
#	boundCRS = "+proj=laea +lat_0=48 +lon_0=9 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m",
#	boundCRS = boundCRS,
#	boundaries = boundaries,
  formulaString = as.formula(paste(obsNames[3],"~1")),
	params =  list(thresh = quantile(observations$value,0.9),block=cellsize),
    outputWhat = list(mean=TRUE, variance=TRUE, excprob = 5.9, cumdistr = 5.9, 
		quantile = .1),
    blockWhat = list(fat=7,fatVar=7,blockMax=TRUE,blockMaxVar = TRUE,blockMin=TRUE),
	class="automap"
)

checkSetup(krigingObject)
krigingObject = preProcess(krigingObject)
krigingObject = estimateParameters(krigingObject)
krigingObject = blockPredict(krigingObject)
krigingObject = postProcess(krigingObject)
predictions = krigingObject$predictions

# Send predictions back to Java. Not sure how to deal with this spatial object though...?
summary(krigingObject$outputTable)
