list.maps <-
function(prefix = "*") {
    raster.list <- execGRASS("g.list", type = "raster", pattern = prefix, intern = TRUE)
    vector.list <- execGRASS("g.list", type = "vector", pattern = prefix, intern = TRUE)
    return(list(raster = raster.list, vector = vector.list))
}
