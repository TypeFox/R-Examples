createSphereSegment <-
function(partOfSP,
                              max.radius=100,  #m
                             name="segment",
                             key.entries = as.numeric(partOfSP@data),
                             scalelist=1,
                             do.sqrt = TRUE,
                             fillColor=rainbow(length(key.entries)),
                             fillOpacity=0.7,
                             map="map",
                             strokeColor="#FDA0C0C",
                             strokeOpacity=1,
                             strokeWeight=1,
                             geodesic="null",
                             clickable=TRUE,
                             zIndex=1) {

 center=coordinates(partOfSP)
fillColor<-as.character(substr(fillColor,1,7))

	obj = as(partOfSP, "SpatialPointsDataFrame")
		z = as.numeric(partOfSP@data)
    	# avoid negative values
    	if (min(key.entries)<0 ){
    	ke<-abs(min(key.entries))+ key.entries+mean(key.entries)
    	}else{ke<-key.entries+min(key.entries)}     # no zeros for radius vecor
    	# creating a vector for subgroupes
    	if(do.sqrt){
    	scale.level<- sqrt(ke/(max(ke)) ) }else{scale.level<-ke/(max(ke))}
	radius.level<-max.radius*scale.level
	 # list of radiuses for createSphereCircle
	radius.vector<-   radius.level
	dfi<-cbind(rep(NA,1+length(radius.vector)))
	dfi[1]=0

	if(max(scale.level)==0)
                          {scale.level=0.1
                           scalelist=0.01}

	dfi[2:(length(radius.vector)+1)]=360/sum(scale.level)*  scale.level


	fi= cbind(rep(NA,2*length(radius.vector)) )
	dfi=as.numeric(dfi)
 for (i in (seq(2,length(fi),by=2))){
 fi[i]=sum(dfi[1:(i/2+1)])
 }
 fi[1]=0
   for (i in (seq(3,length(fi),by=2))){
 fi[i]=fi[i-1]
 }





radius.vector=scalelist*max.radius/6378000
# radius<-	spDistsN1(pts=rbind(c(soil.ll@bbox[1,1],soil.ll@bbox[2,1])), pt=c(soil.ll@bbox[1,2],soil.ll@bbox[2,2]), longlat = T)/6378	## radians
              if (clickable!=FALSE)
                  {clickable='true'}else{ clickable='false' }

              if (geodesic!=FALSE){
              geodesic='true'}else{  geodesic='false'}


		lat1 <- (pi/180)* center[2];
	  lng1 <- (pi/180)* center[1];

  paths<-as.list(rep(NA,length(radius.vector)))

 for(ik in (seq(2,length(fi),by=2))){
     radius<-  radius.vector

    coords<-cbind(rep(NA,11),rep(NA,11))
    coords[1,]  <- c(  center[,1], center[,2])
      j<-2
		for (i in seq(fi[ik-1],fi[ik],length.out=10)) {
			tc <- (pi/180)*i
			y <- asin(sin(lat1)*cos(radius)+cos(lat1)*sin(radius)*cos(tc))
			dlng <- atan2(sin(tc)*sin(radius)*cos(lat1),cos(radius)-sin(lat1)*sin(y))
			x <- ((lng1-dlng+pi) %% (2*pi)) - pi

			coords[j,1] <-c(x*(180/pi))
			coords[j,2] <-c(y*(180/pi))
			j<-j+1
		}



       lonlat<-coords
       pomocno=paste('new google.maps.LatLng(',lonlat[1:((length(lonlat)/2-1)),2],
                    ',',lonlat[1:((length(lonlat)/2-1)),1],'),\n',collapse="")

       paths[[ik/2]]=paste('[',pomocno,'new google.maps.LatLng(',
                 lonlat[length(lonlat)/2,2],',',lonlat[length(lonlat)/2,1],')] ')


                 }

   seg<-""
   for(i in (1:(length(paths)-1))){
    seg<- paste(seg,' new google.maps.Polygon({ \n path:',paths[[i]],', \n','map:',
           map,', \n clickable:',clickable,',\n fillColor: "',
           fillColor[i],'",\n strokeColor: "',strokeColor,'", \n strokeOpacity:',
           strokeOpacity,',\n fillOpacity:',fillOpacity,',\n strokeWeight:',
           strokeWeight,',\n geodesic:',geodesic,',\n zIndex:',zIndex,'}), \n',sep="") }


seg<-paste(seg,' new google.maps.Polygon({ \n path:',paths[[length(paths)]],', \n','map:',
           map,', \n clickable:',clickable,',\n fillColor: "',
           fillColor[length(paths)],'",\n strokeColor: "',strokeColor,'", \n strokeOpacity:',
           strokeOpacity,',\n fillOpacity:',fillOpacity,',\n strokeWeight:',
           strokeWeight,',\n geodesic:',geodesic,',\n zIndex:',zIndex,'})] \n',sep="")


seg<-paste('var ',name, '= [', seg,'\n ',sep="")

return(seg)


}
