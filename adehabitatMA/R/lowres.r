lowres <- function(x, np=2, which.fac=NULL, ...)
{
    ## Verifications
    if (is(x, "SpatialGrid"))
        fullgrid(x) = FALSE
    if (!inherits(x, "SpatialPixelsDataFrame"))
          stop("x should be of class \"SpatialPixelsDataFrame\"")
    pfs <- proj4string(x)
    gr <- gridparameters(x)
    if (nrow(gr) > 2)
        stop("x should be defined in two dimensions")
    if ((gr[1, 2] - gr[2, 2])> get(".adeoptions", envir=.adehabitatMAEnv)$epsilon)
        stop("the cellsize should be the same in x and y directions")
    res <- list()

    for (i in 1:(ncol(slot(x,"data")))) {
        nc <- gr[2, 3]
        nr <- gr[1, 3]
        cs <- gr[2, 2]
        maa<-as.image.SpatialGridDataFrame(x[,i])
        y <- maa$z
        typ <- "numeric"
        if (i%in%which.fac)
            typ="factor"

        ## build a smaller matrix, multiple of np (to avoid "half-pixels")
        y<-y[1:(nr-(((nr/np)-floor(nr/np)))*np),
             1:(nc-(((nc/np)-floor(nc/np)))*np)]
        nr<-nrow(y)
        nc<-ncol(y)

        ## recomputes the levels of the map if it is a factor
        if (typ=="factor") {
            repr<- as.numeric(levels(factor(as.vector(y))))
            y <- as.numeric(as.character(factor(y)))
            y <- matrix(y, nrow=nr, ncol=nc)
        }

        ## Replaces the missing values
        y[is.na(y)]<--9999

        ## the future output
        xs<-matrix(0, nrow=nr/np, ncol=nc/np)

        if (typ == "numeric") {
            ## in case of numeric map: computes the average value for the pixel
            mat<-.C("regrouascnumr", as.double(t(y)), as.double(t(xs)),
                    as.double(nrow(y)), as.double(ncol(y)),
                    as.double(nrow(xs)), as.double(ncol(xs)),
                    PACKAGE = "adehabitatMA")[[2]]
        } else {
            ## in case of factor maps: computes the most frequent value
            ## for the pixel
            mat<-.C("regroufacascr", as.double(t(y)),
                    as.double(t(xs)), as.integer(np),
                    as.integer(length(repr)), as.integer(nrow(y)),
                    as.integer(ncol(y)),
                    as.integer(nrow(xs)), as.integer(ncol(xs)),
                    PACKAGE = "adehabitatMA")[[2]]
        }

        ## The output
        mat<-matrix(mat,ncol=ncol(xs), byrow=TRUE)
        mat[mat==-9999]<-NA
        maa$z <- mat
        maa$x <- mean(maa$x[1:np]) + c(0:(nr/np - 1)) * cs*np
        maa$y <- mean(maa$y[1:np]) + c(0:(nc/np - 1)) * cs*np
        maa <- image2Grid(maa)
        maa <- as(maa, "SpatialPixelsDataFrame")
        gridded(maa) <- TRUE
        res[[i]] <- maa
    }
    names(res) <- names(slot(x, "data"))
    re <- do.call("data.frame",lapply(res, function(x) x[[1]]))
    coordinates(re) <- coordinates(res[[1]])
    gridded(re) <- TRUE
    if (!is.na(pfs))
        proj4string(re) <- CRS(pfs)

    return(re)
}

