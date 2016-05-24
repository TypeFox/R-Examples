### R code from vignette source 'RapidPolygonLookup.Rnw'

###################################################
### code chunk number 1: RapidPolygonLookup.Rnw:69-77
###################################################
library(RgoogleMaps)
library(RapidPolygonLookup)
data(sf.crime.2012)
sf.map <- GetMap(center=c(lat = 37.77173, lon = -122.4306), destfile = "sanFrancisco.png", 
                 size=c(550,550), GRAYSCALE = FALSE, zoom = 13, verbose= FALSE)
color <- ifelse(sf.crime.2012$violent == "TRUE", "red", "green4")
PlotOnStaticMap(MyMap= sf.map, lat= sf.crime.2012$Y, lon= sf.crime.2012$X, 
                col= color, cex= 0.35, pch= 20)


###################################################
### code chunk number 2: RapidPolygonLookup.Rnw:91-93
###################################################
data(california.tract10)
plot(california.tract10)


###################################################
### code chunk number 3: RapidPolygonLookup.Rnw:101-106
###################################################
sf.polys <- CropSpatialPolygonsDataFrame(x= california.tract10, 
                                         bb= data.frame(X=c(-122.5132, -122.37), 
                                                        Y= c(37.70760, 37.81849)))
str(sf.polys, max.level =1)
plotPolys(sf.polys$polys)


###################################################
### code chunk number 4: RapidPolygonLookup.Rnw:116-123
###################################################
XY.kdtree <- RapidPolygonLookup(sf.crime.2012[,c("X","Y")], poly.list= sf.polys, 
                                k= 10, N= 5000, 
                                poly.id= "fips", poly.id.colname= "census.block", 
                                keep.data= TRUE, verbose= FALSE)
table(XY.kdtree$XY$rank, useNA= "always")
hist(XY.kdtree$XY$rank, xlab = "Rank of neighbor", 
     main= "Histogram of number of polygons searched")


