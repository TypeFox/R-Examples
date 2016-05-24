"distfacmap" <- function(x, lev = NULL)
{
    ## Verifications
    if (is(x, "SpatialGrid"))
        fullgrid(x) = FALSE
    if (!inherits(x, "SpatialPixelsDataFrame"))
        stop("x should be of class \"SpatialPixelsDataFrame\"")

    gr <- gridparameters(x)
    pfs <- proj4string(x)
    if (nrow(gr) > 2)
        stop("x should be defined in two dimensions")
    if ((gr[1, 2] - gr[2, 2])> get(".adeoptions", envir=.adehabitatMAEnv)$epsilon)
        stop("the cellsize should be the same in x and y directions")


    ## Bases for the function
    xyc <- as.data.frame(coordinates(x))
    if (is.null(lev)) {
        lev <- levels(factor(x[[1]]))
    } else {
        if (length(levels(factor(x[[1]]))) != length(lev))
            stop("non convenient length for lev")
    }
    li <- list()

    ## For each level of the map:
    for (i in 1:length(lev)) {

        ## keeps only the coordinates of the pixels
        ## corresponding to this level
        tmp <- x[[1]]
        tmp[tmp!=i] <- NA
        tmp[tmp==i] <- 1
        ptsoui <- xyc[!is.na(c(tmp)),]

        ## these objects are passed to a call to the C function "distxyr",
        ## which computes the distance of each pixel to the nearest pixel
        ## for which the level is i
        toto <- .C("distxyr", as.double(t(as.matrix(xyc))),
                   as.double(t(as.matrix(ptsoui))),
                   as.integer(nrow(xyc)), as.integer(nrow(ptsoui)),
                   double(nrow(xyc)), PACKAGE="adehabitatMA")
        li[[i]] <- toto[[5]]
    }

    li <- do.call("data.frame", li)
    names(li) <- paste("level", lev, sep=".")
    coordinates(li) <- xyc
    gridded(li) <- TRUE

    if (!is.na(pfs))
        proj4string(li) <- CRS(pfs)

    return(li)
}

