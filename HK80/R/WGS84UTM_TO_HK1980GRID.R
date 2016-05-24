WGS84UTM_TO_HK1980GRID <-
function(N, E, zone){
    temp1 <- WGS84UTM_TO_WGS84GEO(N, E, zone)
    temp2 <- WGS84GEO_TO_HK80GEO(temp1$latitude, temp1$longitude)
    res <- HK80GEO_TO_HK1980GRID(temp2$latitude, temp2$longitude)
    return(res)
}
