readVista<-function (i="garmin",f = "usb:", type = "w", SPDF=NULL, invisible = TRUE){
	
  if (type=="w"){
        gpsdata <- system(paste("gpsbabel -w -i ",i," -f ",f," -o gpx -F -", sep = ""), intern = TRUE, invisible = TRUE)
				if (!any(grep("<gpx", gpsdata))) stop("Cannot read GPS: check connexion and the device interface argument 'f' (eg. type, port number, etc.)")
        tmp<-tempfile(pattern="readVista")
        if(is.null(SPDF)){
        cat(gpsdata,file=tmp)
			  gpsdata<-data.frame(readOGR(tmp,"waypoints"))[,c(5,20,21,1)]
        names(gpsdata)<-c("ident","long","lat","altitude")
        attributes(gpsdata)$type<-"Waypoints"
			  attributes(gpsdata)$sysTime<-Sys.time()
			  unlink(tmp)
			  return(gpsdata)
        } else {
          cat(gpsdata,file=SPDF)
        }
    }

    if (type=="t"){
        gpsdata <- system(paste("gpsbabel -t -i ",i," -f ",f," -o gpx -F -", sep = ""), intern = TRUE, invisible = TRUE)
				if (!any(grep("<gpx", gpsdata))) stop("Cannot read GPS: check connexion and the device interface argument 'f' (eg. type, port number, etc.)")
        tmp<-tempfile(pattern="readVista")
				if(is.null(SPDF)){
        cat(gpsdata,file=tmp)
				gpsdata<-data.frame(readOGR(tmp,"track_points"))[,c(1,25,26,4)]
				names(gpsdata)<-c("ident","long","lat","altitude")
        gpsdata[,1]<-1:nrow(gpsdata)
				attributes(gpsdata)$type<-"Track"
				attributes(gpsdata)$sysTime<-Sys.time()
				unlink(tmp)
				return(gpsdata)
        } else {
				  cat(gpsdata,file=SPDF)
				}
				
    }
    
}



