HK80GEO_TO_WGS84UTM <-
function(latitude, longitude){
    #### The latitude and longitude should be both in decimal format. 
    temp <- HK80GEO_TO_HK80UTM(latitude, longitude)
    res <- HK80UTM_TO_WGS84UTM(temp$N, temp$E, zone = temp$zone)
    return(res)
}
