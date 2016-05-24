library(rtop)
library(rgdal)
set.seed(1)


# Read data sets
rpath = system.file("extdata",package="rtop")
observations = readOGR(rpath, "observations")
predictionLocations = readOGR(rpath, "predictionLocations")
observations$obs = observations$QSUMMER_OB/observations$AREASQKM

# Create an rtop-object and fit a variogram model
rtopObj = createRtopObject(observations, predictionLocations,
formulaString = obs~1,params = list(gDist = TRUE, rresol = 25))
rtopObj = rtopFitVariogram(rtopObj)

# Exploratory data analyses and visualization of model fit
rtopObj = checkVario(rtopObj, cloud = TRUE, identify = TRUE,
   acor = 0.000001)
rtopObj = checkVario(rtopObj, acor = 0.000001, acomp =
       data.frame(acl1 = c(2,2,2,2,3,3,3,4,4), acl2 = c(2,3,4,5,3,4,5,4,5)))

# Cross-validation
rtopObj = rtopKrige(rtopObj,cv=TRUE)
predictions = rtopObj$predictions
summary(predictions)
sstot = sum((predictions$obs-mean(predictions$obs))^2)
rtopsserr = sum((predictions$obs- predictions$var1.pred)^2)
rtoprsq = 1-rtopsserr/sstot
rtoprsq

# Predictions at new locations - visualize the result on a stream network
rtopObj = rtopKrige(rtopObj)
rnet = readOGR(rpath, "riverNetwork")
pred = rtopObj$predictions
rnet$pred = pred$var1.pred[match(rnet$EZGA, pred$EZGID)]
spplot(rnet, "pred", col.regions = bpy.colors())
at = seq(0,max(rnet$pred,na.rm = TRUE),0.01)
cols = bpy.colors(length(at))
cobs = observations@data[,c("XSTATION", "YSTATION", "obs")]
names(cobs) = c("x","y","obs")
coordinates(cobs) = ~x+y
cobs$class = findInterval(cobs$obs, at)

spplot(rnet,"pred",col.regions = bpy.colors(), at = at, panel = function(x,y, ...){
panel.polygonsplot(x,y, ...)
sp.points(cobs[,"obs"], cex=1, pch = 16, col = cols[cobs$class])
})

#writeOGR(rnet, dsn, layer, "ESRI Shapefile")