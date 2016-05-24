library(automap)
library(psgp)

data(meuse)
observations = data.frame(x = meuse$x,y = meuse$y,value = log(meuse$zinc))
coordinates(observations) = ~x+y
set.seed(13531)
predictionLocations = spsample(observations, 50, "regular")

krigingObject = createIntamapObject(
	observations = observations,
	predictionLocations = predictionLocations,
  formulaString = as.formula(value~1),
	params =  list(doAnisotropy = TRUE, thresh = quantile(observations$value,0.9)),
  outputWhat = list(mean=TRUE, variance=TRUE, excprob = 5.9, cumdistr = 5.9, 
		quantile = .1)
)
class(krigingObject) = c("psgp")

checkSetup(krigingObject)
krigingObject = preProcess(krigingObject)
krigingObject = estimateParameters(krigingObject)
krigingObject = spatialPredict(krigingObject)
krigingObject = postProcess(krigingObject)

# Send predictions back to Java. Not sure how to deal with this spatial object though...?
summary(krigingObject$outputTable)
summary(krigingObject$observations)
summary(autoKrige(value~1,krigingObject$observations,predictionLocations)$krige_output)
autofitVariogram(value~1,krigingObject$observations)$var_model
