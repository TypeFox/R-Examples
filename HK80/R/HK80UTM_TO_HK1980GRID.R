HK80UTM_TO_HK1980GRID <-
function(N, E, zone){
    ### The unit for N and E is meter. The zone is either 49 or 50
    temp <- HK80UTM_TO_HK80GEO(N, E, zone)
    res <- HK80GEO_TO_HK1980GRID(temp$latitude, temp$longitude) 
    return(res)
}
