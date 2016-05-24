createPolygon <-
function(oneSPpolygonsSlot,
                         name="polygon",
                         fillColor="#00AAFF",
                         fillOpacity=0.5,
                         map="map",
                         strokeColor="#FFAA00",
                         strokeOpacity=1,
                         strokeWeight=1,
                         geodesic=TRUE,
                         clickable=TRUE,
                         zIndex="null") {

              if (clickable!=FALSE)
                  {clickable='true'}else{ clickable='false' }

              if (geodesic!=FALSE)
                 {geodesic='true' }else{  geodesic='false'}

listOfPolygonsPerPolygon<-lapply(slot(oneSPpolygonsSlot, "Polygons"),function(x) slot(x, "coords"))
holes<-sapply(slot(oneSPpolygonsSlot, "Polygons"), function(x) slot(x, "ringDir")) 
plotOrder<-slot(oneSPpolygonsSlot, "plotOrder")
pts<-rep("",length(plotOrder))

   if (length(plotOrder)>1) {

              for(j in 1:(length(plotOrder)-1)){
              lonlat<-listOfPolygonsPerPolygon[[plotOrder[j]]]
              #if(holes[j]==-1) {lonlat<-lonlat[ nrow(lonlat):1, ]}
              
              pts[j]=paste('new google.maps.LatLng(',
                           lonlat[1:((length(lonlat)/2-1)),2],',',
                           lonlat[1:((length(lonlat)/2-1)),1],'),\n',collapse="")
              pts[j]=paste(pts[j],'new google.maps.LatLng(',
                           lonlat[length(lonlat)/2,2],',',lonlat[length(lonlat)/2,1],')')
              pts[j]=paste('[',pts[j],'], \n')
              }
       lonlat<-listOfPolygonsPerPolygon[[plotOrder[j+1]]]
       pts[j+1]=paste('new google.maps.LatLng(',
                       lonlat[1:((length(lonlat)/2-1)),2],',',
                       lonlat[1:((length(lonlat)/2-1)),1],'),\n',collapse="")
       pts[j+1]=paste(pts[j+1],'new google.maps.LatLng(',
                     lonlat[length(lonlat)/2,2],',',lonlat[length(lonlat)/2,1],')')
       pts[j+1]=paste('[',pts[j+1],'] \n')
       xxx=paste(pts[1:(j+1)],collapse="")
       paths=paste('[',xxx,']')
         }else{

       lonlat<-slot(slot(oneSPpolygonsSlot, "Polygons")[[1]],"coords")
       paths=paste('new google.maps.LatLng(',
                   lonlat[1:((length(lonlat)/2-1)),2],',',
                   lonlat[1:((length(lonlat)/2-1)),1],'),\n',collapse="")
       paths=paste('[',paths,'new google.maps.LatLng(',
                   lonlat[length(lonlat)/2,2],',',lonlat[length(lonlat)/2,1],')]')

  }






x<-paste('var ',name,'= new google.maps.Polygon({ \n paths:',paths,', \n','map:',
         map,', \n clickable:',clickable,',\n fillColor: "',fillColor,
         '",\n strokeColor: "',strokeColor,'", \n strokeOpacity:',
         strokeOpacity,',\n fillOpacity:',fillOpacity,',\n strokeWeight:',
         strokeWeight,',\n geodesic:',geodesic,',\n zIndex:',zIndex,'});',sep="")

return(x)

}
