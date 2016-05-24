HK80UTM_TO_WGS84UTM <-
function(N, E, zone = c(49, 50)){
    ### The unit for N and E is meter. The zone is either 49 or 50
    if(zone == 49){
        res.N <- N - 195
        res.E <- E + 245
    }
    if(zone == 50){
        res.N <- N - 205
        res.E <- E + 260
    }
    return(list(N = res.N, E = res.E, zone = zone))
}
