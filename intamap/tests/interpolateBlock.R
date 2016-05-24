library(intamap)
data(meuse)
coordinates(meuse) = ~x+y
data(meuse.grid)
coordinates(meuse.grid) = ~x+y
set.seed(13531)

predictionLocations = spsample(meuse,50,"regular")
gridded(predictionLocations) = TRUE
cs = predictionLocations@grid@cellsize[1]/2
meuse$value = log(meuse$zinc)

outputWhat = list(mean=TRUE,variance=TRUE,quantile=0.025,quantile=0.0975)
res1 = interpolateBlock(meuse,predictionLocations,outputWhat,methodName = "automap")$outputTable
summary(t(res1))

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

res2 = interpolateBlock(meuse,predictionLocations,outputWhat,methodName="automap")$outputTable
summary(t(res2))

max((res2-res1)/res1)

