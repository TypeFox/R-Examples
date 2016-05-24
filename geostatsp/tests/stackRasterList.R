library('geostatsp')

grd <- GridTopology(c(1,1), c(1,1), c(10,10))
polys <- as(grd, "SpatialPolygons")
centroids <- coordinates(polys)
x <- centroids[,1]
y <- centroids[,2]
z <- 1.4 + 0.1*x + 0.2*y + 0.002*x*x
xpoly <- SpatialPolygonsDataFrame(polys,
    data=data.frame(x=x, y=y, z=z, row.names=row.names(polys)))

names(xpoly)=paste("stuff", 
    1:length(names(xpoly)), sep="")


template = squareRaster(xpoly, 100)

thebrick = spdfToBrick(
    x=xpoly,
    template=template,
    pattern='^stuff[[:digit:]]$'
    )
    
plot(thebrick[[1]])    
