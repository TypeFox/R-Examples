MapGrid <- function(rast, ...) {
#     layout(matrix(c(1), 1, 1))
#     par(mar = c(4, 4, 4, 4))
    loadNamespace("raster")
    plot(rast, ...)  #, col = colo)
    maps::map("world", add = T, col = "grey")
    
} 
