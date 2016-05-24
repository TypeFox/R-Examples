HK1980GRID_TO_WGS84UTM <-
function(N, E){
    #### The unit for N and E is meter
    temp <- HK1980GRID_TO_HK80GEO(N, E)
    temp2  <- HK80GEO_TO_HK80UTM(temp$latitude, temp$longitude)
    res <- HK80UTM_TO_WGS84UTM(temp2$N, temp2$E, zone = temp2$zone)
    return(res)
}
