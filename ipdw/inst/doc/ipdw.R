### R code from vignette source 'ipdw.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: ipdw00
###################################################
options(prompt = "R> ", continue = "+", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: ipdw-1
###################################################
library("ipdw")
library("geoR")
data(kattegat)
katproj<-c("+proj=utm +zone=32 +ellps=GRS80 +units=m +no_defs")

pols1<-Polygons(list(Polygon(kattegat$dk[1])),"pol1")
pols2<-Polygons(list(Polygon(kattegat$dk[2])),"pol2")
pols3<-Polygons(list(Polygon(kattegat$dk[3])),"pol3")
pols4<-Polygons(list(Polygon(kattegat$dk[4])),"pol4")
pols5<-Polygons(list(Polygon(kattegat$dk[5])),"pol5")
pols6<-Polygons(list(Polygon(kattegat$dk[6])),"pol6")
pols7<-Polygons(list(Polygon(kattegat$dk[7])),"pol7")
pols8<-Polygons(list(Polygon(kattegat$dk[8])),"pol8")
pols9<-Polygons(list(Polygon(kattegat$dk[9])),"pol9")
pols10<-Polygons(list(Polygon(kattegat$dk[10])),"pol10")
pols11<-Polygons(list(Polygon(kattegat$dk[11])),"pol11")
pols12<-Polygons(list(Polygon(kattegat$dk[12])),"pol12")
pols<-SpatialPolygons(list(pols1,pols2,pols3,pols4,pols5,pols6,
                           pols7,pols8,pols9,pols10,pols11,pols12),1:12)


###################################################
### code chunk number 3: ipdw-2
###################################################
projection(pols)<-katproj
costras<-costrasterGen(kattegat$coords,pols,extent="pnts",katproj)
#insert contiguous barrier
costras[160:170,1:80] <- 10000


###################################################
### code chunk number 4: ipdw-3
###################################################
#find average nearest neighbor
library(spatstat)
W=owin(range(kattegat$coords[,1]),range(kattegat$coords[,2]))
kat.pp<-ppp(kattegat$coords[,1],kattegat$coords[,2],window=W)
mean.neighdist<-mean(nndist(kat.pp))

#grid building
gridsize<-mean.neighdist*2
grainscale.fac<-gridsize/res(costras)[1]
gridras<-aggregate(costras,fact=grainscale.fac)
gridpol<-rasterToPolygons(gridras)
gridpol$value<-row.names(gridpol)

#spatial join
kat.df<-data.frame(kattegat)
coordinates(kat.df)<-~x.utm+y.utm
projection(kat.df)<-katproj
fulldataset.over<-over(kat.df,gridpol)
fulldataset.over<-cbind(data.frame(fulldataset.over),data.frame(kat.df))

#grid selection
set.seed(2)
gridlev<-unique(fulldataset.over$value)
for(i in 1:length(gridlev)){
  activesub<-subset(fulldataset.over,fulldataset.over$value==gridlev[i])
  selectnum<-gdata::resample(1:nrow(activesub),1)
  if(i==1){
    training<-activesub[selectnum,]
  }
  else{
    training<-rbind(training,activesub[selectnum,])
  }
}


###################################################
### code chunk number 5: ipdw-4
###################################################
validate<-fulldataset.over[!(row.names(fulldataset.over) %in% row.names(training)),]
xy<-cbind(training$x.utm,training$y.utm)
training<-SpatialPointsDataFrame(xy,training)
xy<-cbind(validate$x.utm,validate$y.utm)
validate<-SpatialPointsDataFrame(xy,validate)
projection(training)<-katproj
projection(validate)<-katproj


###################################################
### code chunk number 6: figure0
###################################################
plot(costras)
points(training)
points(validate,col="red")


###################################################
### code chunk number 7: ipdw-5
###################################################
paramlist <- c("data")
final.ipdw <- ipdw(training, costras, range = mean.neighdist * 10, paramlist,
									 overlapped = TRUE)


###################################################
### code chunk number 8: figure1
###################################################
plot(final.ipdw, main = "Kattegat salinity (ppt)")


###################################################
### code chunk number 9: ipdw-6
###################################################
idw.grid<-rasterToPoints(costras,fun=function(x){x<10000},spatial=TRUE)
gridded(idw.grid)=TRUE
kat.idw<-gstat::idw(data~1,training,idw.grid,maxdist=mean.neighdist*10,debug.level=0)
final.idw<-raster(kat.idw)


###################################################
### code chunk number 10: figure2
###################################################
par(mfrow=c(1,3))
plot(final.ipdw,main="IPDW")
plot(final.idw,main="IDW")
plot(final.idw-final.ipdw,main="IDW versus IPDW")


###################################################
### code chunk number 11: ipdw-7
###################################################
measured.spdf<-data.frame(validate$data)
coordinates(measured.spdf)<-coordinates(validate)

valid.ipdw<-errorGen(final.ipdw,measured.spdf,measured.spdf@data)
valid.idw<-errorGen(final.idw,measured.spdf,measured.spdf@data)


###################################################
### code chunk number 12: figure3
###################################################
par(mfrow=c(1,2))
valid.ipdw<-errorGen(final.ipdw,measured.spdf,measured.spdf@data,plot=TRUE,title="IPDW")
valid.idw<-errorGen(final.idw,measured.spdf,measured.spdf@data,plot=TRUE,title="IDW")


