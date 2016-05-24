WGS84UTM_TO_HK80GEO <-
function(N, E, zone){
    ### The unit for N and E is meter, and zone is either 49 or 50
    temp1 <- WGS84UTM_TO_WGS84GEO(N, E, zone)
    res <- WGS84GEO_TO_HK80GEO(temp1$latitude, temp1$longitude)
    return(res)
}
