
uploadGPS<-function(gpx, f="usb:",type="w") {
  if (type=="w" | type=="t"){
  if (Sys.info()["sysname"]!="Windows") {# OK sous linux
        system(paste("gpsbabel -",type," -i gpx -f ",gpx," -o garmin -F ",f,sep=""))
        }
    else {# OK sous Windows
         system(paste("gpsbabel -",type," -i gpx -f ",gpx," -o garmin -F ",f,sep=""))
    }
  } else stop("type must be w (waypoints) or t (track)")
}