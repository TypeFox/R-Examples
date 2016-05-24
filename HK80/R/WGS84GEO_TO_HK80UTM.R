WGS84GEO_TO_HK80UTM <-
function(latitude, longitude){
    #### The latitude and longitude should be both in decimal format. 
    temp <- WGS84GEO_TO_WGS84UTM(latitude, longitude)
    res <- WGS84UTM_TO_HK80UTM(temp$N, temp$E, temp$zone)
    return(res)
}
