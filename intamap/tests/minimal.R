library(intamap)

# set up data:
data(meuse)
coordinates(meuse) = ~x+y
meuse$value = log(meuse$zinc)
data(meuse.grid)
gridded(meuse.grid) = ~x+y

# set up intamap object:
obj = createIntamapObject(
	observations = meuse,
	predictionLocations = meuse.grid,
	constantBias = 0,
	outputWhat = list(mean = 1, 
		variance = 1, 
		quantile = 0.05, 
		quantile = 0.5, 
		quantile = 0.95,
		excprob = 5.5,
		excprob = 6.6,
		cumdistr = 5.5,
		cumdistr = 6.6,
		cumdistr = 7.9
		)
)
class(obj) = "linearVariogram"

# check:
checkSetup(obj)

# do interpolation steps:
obj = preProcess(obj)
obj = estimateParameters(obj)
obj = spatialPredict(obj)
obj = postProcess(obj)
summary(obj$predictions)
output = obj$outputTable
