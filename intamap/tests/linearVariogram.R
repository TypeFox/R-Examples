library(intamap)

# set up data:
data(meuse)
coordinates(meuse) = ~x+y
meuse$value = log(meuse$zinc)
data(meuse.grid)
gridded(meuse.grid) = ~x+y
proj4string(meuse) = CRS("+init=epsg:28992")
proj4string(meuse.grid) = CRS("+init=epsg:28992")
set.seed(13531)

mgrid = coarsenGrid(meuse.grid,4)
# set up intamap object:
obj = createIntamapObject(
	observations = meuse,
	predictionLocations = mgrid,
	targetCRS = "+init=epsg:3035",
	params = list(predictType=list(quantiles=c(0.05,0.5,0.95)),thresh=c(5.5,6.6)),
  outputWhat = list(mean = 1,
         variance = 1,
          quantile = 0.05,
          quantile = 0.5,
          quantile = 0.95,
      		excprob = 5.5,
      		excprob = 6.6,
      		cumdistr = 5.5,
      		cumdistr = 6.6,
      		cumdistr = 7.9)
)
class(obj) = "linearVariogram"

# check:
checkSetup(obj)

# do interpolation steps:
obj = preProcess(obj)
obj = estimateParameters(obj)
obj = spatialPredict(obj)
obj = postProcess(obj)
output = obj$predictions

# generate some output:
summary(obj$predictions)
gridded(output) = FALSE
summary(output)
