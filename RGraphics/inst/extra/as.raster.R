as.raster.pixmapGrey <- function(x) {
    r <- gray(t(x@grey))
    dim(r) <- x@size
    class(r) <- "raster"
    r
}
as.raster.pixmapRGB <- function(x) {
    r <- rgb(t(x@red), t(x@green), t(x@blue))
    dim(r) <- x@size
    class(r) <- "raster"
    r
}
