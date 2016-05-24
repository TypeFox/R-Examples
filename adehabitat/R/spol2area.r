"spol2area" <- function(sr)
{
    ## Verifications
    if (inherits(sr, "SpatialPolygonsDataFrame"))
        sr <- polygons(sr)
    if (!inherits(sr, "SpatialPolygons"))
        stop("sr should be of class \"SpatialPolygons\" or \"SpatialPolygonsDataFrame\"")

    ## Gets the polygons in the object sr
    pol <- getSpPpolygonsSlot(sr)
    warh <- 0
    warh2 <- 0
    warz <- 0

    ## For each SpatialPolygons
    res <- lapply(pol, function(x) {

        ## gets the polygons and the ID
        y <- getPolygonsPolygonsSlot(x)
        nom <- getPolygonsIDSlot(x)
        ll <- length(y)

        ## Identify the holes
        hh <- unlist(lapply(y, function(o) getPolygonHoleSlot(o)))
        hol <- sum(hh)
        ll <- ll-hol

        ## One polygon
        if (ll == 1) {

            if (hol == 0) {

                ## No hole -> creates the object area
                re <- as.data.frame(getPolygonCoordsSlot(y[[1]]))
                re <- data.frame( fac = factor(rep(nom,length(re[,1]))), re)
                names(re) <- c("fac", "x", "y")
            }

            if (hol != 0) {

                ## Hole present: prevision of the warning
                ## we will delete the hole
                warh <- warh+hol
                warh2 <- warh2+1

                ## Creation of the object area
                re <- as.data.frame(getPolygonCoordsSlot(y[!hh][[1]]))
                re <- data.frame( fac = factor(rep(nom,length(re[,1]))), re)
                names(re) <- c("fac", "x", "y")
            }
        }

        ## More than one polygon
        if (ll > 1) {
            warz <- warz+1

            ## No hole
            if (hol == 0) {

                ## pools the polygoons into the same area object
                nom <- paste(nom, 1:ll, sep=".")

                ## and creates the object of class "area"
                re1 <- lapply(y, function(o) as.data.frame(getPolygonCoordsSlot(o)))
                re <- do.call("rbind.data.frame",
                              lapply(1:length(re1), function(i) {
                                  u <- data.frame(fac=factor(rep(nom[i],
                                                  length(re1[[i]][,1]))),
                                                  re1[[i]])
                                  names(u) <- c("fac", "x", "y")
                                  return(u)
                              }))
            }

            ## Hole present
            if (hol!=0) {

                ## predict the warning
                warh <- warh+hol
                warh2 <- warh2+1

                ## creates the object of class "area"
                nom <- paste(nom, 1:ll, sep=".")
                y <- y[!hh]
                re1 <- lapply(y, function(o) as.data.frame(getPolygonCoordsSlot(o)))
                re <- do.call("rbind.data.frame",
                              lapply(1:length(re1), function(i) {
                                  u <- data.frame(fac=factor(rep(nom[i],
                                                  length(re1[[i]][,1]))),
                                                  re1[[i]])
                                  names(u) <- c("fac", "x", "y")
                                  return(u)
                              }))
            }
        }

        ## Output
        return(list(re,warh2, warh, warz))
    })

    ## Output object of class "area"
    warh2 <- sum(unlist(lapply(res, function(x) x[[2]])))
    warh <- sum(unlist(lapply(res, function(x) x[[3]])))
    warz <- sum(unlist(lapply(res, function(x) x[[4]])))
    res <- lapply(res, function(x) x[[1]])
    res <- do.call("rbind.data.frame", res)
    res <- as.area(res)

    ## and a warning if holes are present
    if (warh2>0) {
        warning(paste("Area objects do not take into account holes in polygon.\n",
                      warh, "holes have been deleted from the data, belonging to\n",
                      warh2, "polygons"))
    }

    ## Output
    return(res)
}

