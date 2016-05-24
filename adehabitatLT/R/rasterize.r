rasterize.ltraj <- function (ltr, map)
{
    if (!inherits(ltr, "ltraj"))
        stop("ltr should be of class ltraj")
    if (!inherits(map, "SpatialPixels"))
        stop("map should inherit the class SpatialPixels")

    pa <- gridparameters(map)
    pfs <- proj4string(map)
    xll <- pa[1, 1]
    yll <- pa[2, 1]
    cs <- pa[1, 2]
    lapply(1:length(ltr), function(i) {
        x <- ltr[[i]][, c("x", "y")]
        index <- NULL
        index2 <- c(1:nrow(x))
        ## which step indices
        for (j in 1:(nrow(x)-1)) {
            if ((!is.na(sum(unlist(x[j,]))))&(!is.na(sum(unlist(x[j+1,])))))
                index <- c(index,j)
        }
        index2 <- index2[!apply(x,1,function(y) any(is.na(y)))]
        x <- x[!apply(x,1,function(y) any(is.na(y))),]
        res <- .Call("RasterPas", x, xll, yll, cs, as.integer(0),
                     PACKAGE = "adehabitatLT")
        res <- as.data.frame(res)
        names(res) <- c("x", "y", "step")
        res$step <- index2[res$step]
        res <- res[res$step%in%index,]
        coordinates(res) <- c("x", "y")
        if (!is.na(pfs))
            proj4string(res) <- CRS(pfs)
        return(res)
    })
}
