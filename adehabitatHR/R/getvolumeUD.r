"getvolumeUD" <- function(x, standardize=FALSE)
{
    ## Verifications
    if ((!inherits(x, "estUDm"))&(!inherits(x, "estUD")))
        stop("x should be an object of class \"estUD\" or \"estUDm\"")
    if (inherits(x, "estUDm")) {
        if (slot(x[[1]], "vol"))
            stop("already a volume under UD")
    } else {
        if (slot(x, "vol"))
            stop("already a volume under UD")
    }

    if (inherits(x, "estUD")) {
        pfs <- proj4string(x)
        sg <- as(x, "SpatialPixelsDataFrame")
        sg <- as(sg, "SpatialGridDataFrame")
        gr <- gridparameters(sg)
        uu <- names(sg)
        gri <- as.image.SpatialGridDataFrame(sg)
        xyok <- expand.grid(gri$y,gri$x)[,2:1]
        asc <- gri$z

        cs <- gr[1, 2]
        if (standardize) {
            asc <- asc/(sum(asc)*cs*cs)
        }

        ## computes the volume for each pixel
        ## thanks to a call to the C function calcvolume
        v<-.C("calcvolume", as.double(t(asc)), as.integer(ncol(asc)),
              as.integer(nrow(asc)), as.double(cs), PACKAGE="adehabitatHR")[[1]]

        ## standardize it so that the total volume is 1 over the area
        index<-1:length(v)
        vord<-v[order(v, decreasing=TRUE)]
        indord<-index[order(v, decreasing=TRUE)]
        vsu<-cumsum(vord)

        ## output
        vreord<-data.frame(n=vsu[order(indord)]*100)
        coordinates(vreord) <- xyok
        gridded(vreord) <- TRUE
        vreord <- new("estUD", vreord)
        slot(vreord, "h") <- x@h
        slot(vreord, "vol") <- TRUE
        proj4string(vreord) <- CRS(pfs)
        return(vreord)

    } else {

        lire <- list()
        for (i in 1:length(x)) {
            lire[[i]] <- getvolumeUD(x[[i]], standardize=standardize)
        }
        names(lire) <- names(x)
        class(lire) <- "estUDm"
        return(lire)

    }
}


image.estUDm <- function(x, ...)
{
    if (!inherits(x, "estUDm"))
        stop("x should be of class \"estUDm\"")
    opar <- par(mfrow=n2mfrow(length(x)))
    tmp <- lapply(1:length(x), function(i) {
        upar <- par(mar=c(0,0,2,0))
        image(x[[i]], ...)
        title(main=names(x)[i])
        box()
        par(upar)
    })
    par(opar)
}
