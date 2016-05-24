gps2gpx<-function (filename = "", i = "garmin", f = "usb:", type = "w", 
    invisible = TRUE) 
{
    if (toupper(substr(filename, nchar(filename) - 3, nchar(filename))) != 
        ".GPX" & filename != "") 
        filename <- paste(filename, ".gpx", sep = "")
    gpsdata <- system(paste("gpsbabel -",type," -i ", i, " -f ", f, 
        " -o gpx -F -", sep = ""), intern = TRUE, invisible = TRUE)
    if (any(grep("Can't init", gpsdata))) 
        stop("Cannot read GPS: check connexion and the device interface argument 'f' (eg. type, port number, etc.)")
    if (filename == "") 
        gpsdata
    else cat(gpsdata, file = filename)
}
