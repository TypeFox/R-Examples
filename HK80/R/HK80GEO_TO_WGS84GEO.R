HK80GEO_TO_WGS84GEO <-
function(latitude, longitude){
    #### The latitude and longitude should be both in decimal format. 
    lat <- latitude - 5.5/3600
    long <- longitude + 8.8/3600
    return(list(latitude = lat, longitude = long))
}
