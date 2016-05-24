HK80UTM_TO_WGS84GEO <-
function(N, E, zone){
    ### The unit for N and E is meter. The zone is either 49 or 50
    temp <- HK80UTM_TO_WGS84UTM(N, E, zone)
    res  <- WGS84UTM_TO_WGS84GEO(temp$N, temp$E, temp$zone)
    return(res)
}
