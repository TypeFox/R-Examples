createSphereCircle <-
function(center,
                             radius=20,  #m
                             name="polygon",
                             fillColor="#00AAFF",
                             fillOpacity=0.7,
                             map="map",
                             strokeColor="#FDA0C0C",
                             strokeOpacity=1,
                             strokeWeight=1,
                             geodesic="null",
                             clickable=TRUE,
                             zIndex=1) {

radius=radius/6378000
# radius<-	spDistsN1(pts=rbind(c(soil.ll@bbox[1,1],soil.ll@bbox[2,1])), pt=c(soil.ll@bbox[1,2],soil.ll@bbox[2,2]), longlat = T)/6378	## radians

              if (clickable!=FALSE)
                  clickable='true'
              else{
                  clickable='false' }


              if (geodesic!=FALSE)
                 geodesic='true'
              else{
                 geodesic='false'}


		lat1 <- (pi/180)* center[2];
	  lng1 <- (pi/180)* center[1];

    coords<-cbind(rep(NA,25),rep(NA,25))
      j<-1
		for (i in seq(0,360,by=15)) {
			tc <- (pi/180)*i
			y <- asin(sin(lat1)*cos(radius)+cos(lat1)*sin(radius)*cos(tc))
			dlng <- atan2(sin(tc)*sin(radius)*cos(lat1),cos(radius)-sin(lat1)*sin(y))
			x <- ((lng1-dlng+pi) %% (2*pi)) - pi

			coords[j,1] <-c(x*(180/pi))
			coords[j,2] <-c(y*(180/pi))
			j<-j+1
		}



       lonlat<-coords
       paths=paste('new google.maps.LatLng(',lonlat[1:((length(lonlat)/2-1)),2],
                      ',',lonlat[1:((length(lonlat)/2-1)),1],'),\n',collapse="")
       paths=paste('[',paths,'new google.maps.LatLng(',
                  lonlat[length(lonlat)/2,2],',',lonlat[length(lonlat)/2,1],')]')

 x<-paste('var ',name,'= new google.maps.Polygon({ \n path:',paths,
          ', \n','map:',map,', \n clickable:',clickable,',\n fillColor: "',
          fillColor,'",\n strokeColor: "',strokeColor,'", \n strokeOpacity:',
          strokeOpacity,',\n fillOpacity:',fillOpacity,',\n strokeWeight:',
          strokeWeight,',\n geodesic:',geodesic,',\n zIndex:',zIndex,'});',sep="")

return(x)


}
