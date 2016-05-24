
QGIS2sp<-function(df=FALSE){

db<-read.delim("clipboard")

mylst<-rep(list(NA),nrow(db))
for(i in 1:nrow(db)) {
  mylst[[i]]<-readWKT(db[i,1])
  if (class(mylst[[i]])!="SpatialPoints") mylst[[i]]<-spChFIDs(mylst[[i]],as.character(i))
} 

mytop<-mylst[[1]]
for(i in 2:length(mylst)) mytop<-spRbind(mytop,mylst[[i]])

if(!df) {
  type<-class(mytop)
  switch(type,
       SpatialPoints=SpatialPointsDataFrame(mytop,db[,2:length(db),drop=FALSE]),
       SpatialLines=SpatialLinesDataFrame(mytop,db[,2:length(db),drop=FALSE]),
       SpatialPolygons=SpatialPolygonsDataFrame(mytop,db[,2:length(db),drop=FALSE])
  )
  } else {
    dfr<-data.frame(coordinates(mytop),db[,2:length(db),drop=FALSE],row.names = 1:nrow(db))
      names(dfr)[1:2]<-c("long","lat")
      return(dfr)
  }
}