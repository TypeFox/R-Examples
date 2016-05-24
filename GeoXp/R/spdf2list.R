`spdf2list` <-
function(data)
{
if(class(data)[1]!="SpatialPolygonsDataFrame")
stop("Must be a SpatialPolygonsDataFrame object")

poly<-data@polygons
n<-length(poly)

X<- poly[1][[1]]@Polygons[[1]]@labpt[1]
Y<- poly[1][[1]]@Polygons[[1]]@labpt[2]
contours<-rbind(NA,NA,NA,poly[1][[1]]@Polygons[[1]]@coords)

 for(i in 2:n)
  {  
   X<-rbind(X,poly[i][[1]]@Polygons[[1]]@labpt[1])
   Y<-rbind(Y,poly[i][[1]]@Polygons[[1]]@labpt[2])
   for(k in 1:length(poly[i][[1]]@Polygons))
    {contours=rbind(contours,NA,NA,NA,poly[i][[1]]@Polygons[[k]]@coords)
    }
  }
contours=rbind(contours,NA,NA,NA)

return(list(X=X,Y=Y,poly=contours))

}

