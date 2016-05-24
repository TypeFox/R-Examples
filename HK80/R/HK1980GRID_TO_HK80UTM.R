HK1980GRID_TO_HK80UTM <-
function(N, E){
     ### The unit for N and E is meter
     temp <- HK1980GRID_TO_HK80GEO(N, E)
     res <- HK80GEO_TO_HK80UTM(temp$latitude, temp$longitude)
     return(res)
}
