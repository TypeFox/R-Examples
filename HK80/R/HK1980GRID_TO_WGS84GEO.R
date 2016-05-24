HK1980GRID_TO_WGS84GEO <-
function(N, E){
    ### The unit for N and E is meters
    temp <- HK1980GRID_TO_HK80GEO(N, E)
    res  <- HK80GEO_TO_WGS84GEO(temp$latitude, temp$longitude)
    return(res)
}
