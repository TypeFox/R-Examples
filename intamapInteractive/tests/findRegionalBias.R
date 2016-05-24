set.seed(1)
library(intamapInteractive)

data(meuse)
observations = data.frame(x = meuse$x,y = meuse$y,value = log(meuse$zinc))
coordinates(observations) = ~x+y
pBoundaries = spsample(observations, 10, "regular",bb = bbox(observations) +  
              matrix(c(-400,-400,400,400),ncol=2),offset=c(0,0))
gridded(pBoundaries) = TRUE
cs = pBoundaries@grid@cellsize[1]/2

Srl = list()
nb = dim(coordinates(pBoundaries))[1]
for (i in 1:nb) {
  pt1 = coordinates(pBoundaries)[i,]
  x1 = pt1[1]-cs
  x2 = pt1[1]+cs
  y1 = pt1[2]-cs
  y2 = pt1[2]+cs

  boun = data.frame(x=c(x1,x2,x2,x1,x1),y=c(y1,y1,y2,y2,y1))
  coordinates(boun) = ~x+y
  boun = Polygon(boun)
  Srl[[i]] = Polygons(list(boun),ID = as.character(i))
}
pBoundaries = SpatialPolygonsDataFrame(SpatialPolygons(Srl),
                                      data = data.frame(ID=c(1:nb)))
observations$ID = over(observations, geometry(pBoundaries))
blines = findBoundaryLines(pBoundaries, regCode = "ID")
rb = findRegionalBias(observations, blines, value~1, regCode = "ID")
rb$regionalBias

obs2 = removeRegionalBias(observations, rb, value~1, regCode = "ID")

summary(obs2)
