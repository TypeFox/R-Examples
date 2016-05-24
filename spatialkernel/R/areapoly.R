## area of a polygon
## "." in name string may be confused with class method???
areapoly <- function(poly) {
    npoly <- nrow(poly)
    asign <- 1
    ans <- .C("area_poly", as.double(poly), as.integer(npoly),
              area=as.double(0), PACKAGE="spatialkernel")$area
    if(ans < 0) {
        asign <- -1
        ans <- -ans
    }
    invisible(list(area=ans, sign=asign, poly=poly))
}
