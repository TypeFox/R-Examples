kver2spol <- function(kv)
{
    x <- kv
    if (!inherits(x, "kver"))
        stop("x should be of class \"kver\"")
    lipols <- lapply(1:length(x), function(i) {
        y <- x[[i]]
        class(y) <- c("data.frame","list")
        res <- split(y[,2:3],y[,1])
        lipol <- lapply(res, function(z) {
            if (sum(abs(z[1,] - z[nrow(z),]))>1e-16)
                z <- rbind(z,z[1,])
            Polygon(as.matrix(z))
        })
        pols <- Polygons(lipol, ID = names(x)[i])
        return(pols)
    })
    return(SpatialPolygons(lipols))
}



