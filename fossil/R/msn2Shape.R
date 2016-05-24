`msn2Shape` <-
function(msn,lats,dist=NULL) {
    require(shapefiles)
    #tests if lats are a matrix or spatial points
    if (class(lats)=='SpatialPoints') lats<-coordinates(lats)
    dis <- is.null(dist)
    n <- dim(msn)[1]
    nm <- matrix(0, n, n)
    nm[lower.tri(nm)] <- msn[lower.tri(msn)]
    l <- 2 * sum(nm)
    if (dis == TRUE) 
        b <- 3
    else b <- 4
    a <- matrix(0, l, 3)
    a[, 1] <- rep(1:(l/2), each = 2)
    if (dis == TRUE) 
        aa <- data.frame(Id = 1:(l/2), localityA = NA, localityB = NA)
    else aa <- data.frame(Id = 1:(l/2), localityA = NA, localityB = NA, 
        dist = NA)
    colnames(a) <- c("Id", "X", "Y")
    rownames(a) <- 1:l
    ind <- colnames(msn)
    li <- 1
    for (i in 2:n) {
        for (j in 1:(i - 1)) {
            if (nm[i, j] == 1) {
                a[li, 2] <- lats[ind[i], 1]
                a[li, 3] <- lats[ind[i], 2]
                li <- li + 1
                a[li, 2] <- lats[ind[j], 1]
                a[li, 3] <- lats[ind[j], 2]
                aa[(li/2), 2] <- ind[i]
                aa[(li/2), 3] <- ind[j]
                if (dis == FALSE) 
                  aa[(li/2), 4] <- as.matrix(dist)[i, j]
                li <- li + 1
            }
        }
    }
    a <- as.data.frame(a)
    ashp <- convert.to.shapefile(a, aa, "Id", 3)
    return(ashp)
}
