WGS84GEO_TO_HK1980GRID <-
function(latitude, longitude){
    #### The latitude and longitude should be in decimal format. 
    temp <- WGS84GEO_TO_HK80GEO(latitude, longitude)
    res <- HK80GEO_TO_HK1980GRID(temp$latitude, temp$longitude)
    return(res)
}
