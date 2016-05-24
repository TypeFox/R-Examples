elSPDF <-function(SPDF,
                  zcol=1:3,
                  scale_e=10){
pointsdf=SPDF
pointsdf<-pointsdf[which(pointsdf@data[,zcol[1]]!=0 & pointsdf@data[,zcol[2]]!=0),]

Pols=lapply(1:length(pointsdf@data[,1]), function (i){
ellip(a=pointsdf@data[i,zcol[1] ]*scale_e,b=pointsdf@data[i,zcol[2] ]*scale_e, alpha = pointsdf@data[i,zcol[3] ], 
         loc = c(pointsdf@coords[i,1],pointsdf@coords[i,2]), n =100,Id=paste(i) ) } )
            
SPl<-SpatialPolygons(Pols,proj4string =pointsdf@proj4string)
SPldf<-SpatialPolygonsDataFrame(SPl,pointsdf@data[which(pointsdf@data[,zcol[1]]!=0 & pointsdf@data[,zcol[2]]!=0),],match.ID = FALSE)
return(SPldf)
}
