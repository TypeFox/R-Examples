`lats2Shape` <-
function(lats) {
  require(shapefiles)
  #tests if lats are a matrix or spatial points
  if (class(lats)=='SpatialPoints') lats<-coordinates(lats)
  n<-dim(lats)[1]
  a<-data.frame(Id=1:n,X=lats[,1],Y=lats[,2])    
  aa<-data.frame(Id=1:n,locality=row.names(as.data.frame(lats)),lats)
  ashp<-convert.to.shapefile(a,aa,"Id",1)
  return(ashp)
}

