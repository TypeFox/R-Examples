"area2spol" <- function(ar)
{
    ## Verifications
    if (!inherits(ar, "area"))
        stop("ar should be of class \"area\"")

    ## splits ar into a list where each element is a polygon
    class(ar) <- "data.frame"
    li <- split(ar[,2:3],ar[,1])

    ## stores the elements as SpatialPolygons
    res <- lapply(li, function(x) {

        ## Verification that the polygon is closed
        if (!all(unlist(x[1,]==x[nrow(x),])))
            x <- rbind(x,x[1,])

        ## converts as spol
        x <- as.matrix(x)
        y <- Polygon(x, hole=FALSE)
        if (y@ringDir<0)
            y <- Polygon(x[nrow(x):1,], hole=FALSE)
        return(y)
    })
    ## The output
    resb <- SpatialPolygons(lapply(1:length(res),
                               function(i) Polygons(list(res[[i]]),
                                                    names(res)[i])))
    return(resb)
  }

