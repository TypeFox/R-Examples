WGS84GEO_TO_HK80GEO <-
function(latitude, longitude){
    #### The latitude and longitude should be both in decimal format. 
    lat <- latitude + 5.5/3600
    long <- longitude - 8.8/3600
    return(list(latitude = lat, longitude = long))
}
