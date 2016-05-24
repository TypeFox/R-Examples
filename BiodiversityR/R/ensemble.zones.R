`ensemble.zones` <- function(
    presence.raster=NULL, centroid.object=NULL, x=NULL, ext=NULL,
    RASTER.species.name=centroid.object$name, RASTER.stack.name = x@title,
    RASTER.format="raster", RASTER.datatype="INT1S", RASTER.NAflag=-127,
    KML.out=FALSE, KML.maxpixels=100000, KML.blur=10
)
{
    .BiodiversityR <- new.env()
#    if (! require(dismo)) {stop("Please install the dismo package")}
    if (is.null(presence.raster) == T) {stop("value for parameter presence.raster is missing (RasterLayer object)")}
    if(inherits(presence.raster, "RasterLayer") == F) {stop("x is not a RasterLayer object")}
    if(is.null(x) == T) {stop("value for parameter x is missing (RasterStack object)")}
    if(inherits(x, "RasterStack") == F) {stop("x is not a RasterStack object")}
    if (is.null(centroid.object) == T) {stop("value for parameter centroid.object is missing (hint: use the ensemble.centroids function)")}

# 
predict.zone <- function(object=centroid.object, newdata=newdata) {
    centroids <- object$centroids
    cov.mahal <- object$cov.mahal
    nc <- nrow(centroids)
    result <- data.frame(array(0, dim=c(nrow(newdata), nc)))
    for (i in 1:nc) {
        result[,i] <- mahalanobis(newdata, center=as.numeric(centroids[i,]), cov=cov.mahal)
    }
    p <- apply(result[, 1:nc], 1, which.min)
    p <- as.numeric(p)
    return(p)
}

#
# check if all variables are present
    vars <- names(centroid.object$centroids)
    vars.x <- names(x)
    nv <- length(vars) 
    for (i in 1:nv) {
        if (any(vars.x==vars[i]) == F) {stop("explanatory variable '", vars[i], "' not among grid layers of RasterStack x \n", sep = "")}
    }
    nv <- length(vars.x) 
    for (i in 1:nv) {
        if (any(vars==vars.x[i]) == F) {
            cat(paste("\n", "NOTE: RasterStack layer '", vars.x[i], "' was not documented in the centroids data set", "\n", sep = ""))
            x <- raster::dropLayer(x, which(names(x) %in% c(vars.x[i]) ))
        }
    }

# same extent for predictors and presence map
    if (is.null(ext) == F) {
        if(length(x@title) == 0) {x@title <- "stack1"}
        title.old <- x@title
        x <- raster::crop(x, y=ext, snap="in")
        x@title <- title.old
        presence.raster <- raster::crop(presence.raster, y=ext, snap="in")
    }

# avoid problems with non-existing directories and prepare for output
    dir.create("ensembles", showWarnings = F)
    dir.create("ensembles/zones", showWarnings = F)
    if(KML.out == T) {
        dir.create("kml", showWarnings = F)
        dir.create("kml/zones", showWarnings = F)
    }
    stack.title <- RASTER.stack.name
    if (gsub(".", "_", stack.title, fixed=T) != stack.title) {cat(paste("\n", "WARNING: title of stack (", stack.title, ") contains '.'", "\n\n", sep = ""))}
    rasterfull <- paste("ensembles/zones/", RASTER.species.name, "_", stack.title , sep="")
    kmlfull <- paste("kml/zones/", RASTER.species.name, "_", stack.title , sep="")

#
# predict
    tryCatch(zones.raster <- raster::predict(object=x, model=centroid.object, fun=predict.zone, na.rm=TRUE, 
            filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format),
        error= function(err) {print(paste("prediction of zones failed"))},
        silent=F)

# mask the presence area, including areas that are NA in presence raster
    zones.raster <- raster::mask(zones.raster, presence.raster, inverse=T, maskvalue=1)
    zones.raster <- raster::mask(zones.raster, presence.raster, inverse=F)
    cat(paste("\n", "raster layer with zones created", "\n", sep = ""))
    print(raster::freq(zones.raster))

#
# avoid possible problems with saving of names of the raster layers
    raster::writeRaster(zones.raster, filename="working.grd", overwrite=T)
    working.raster <- raster::raster("working.grd")
    names(working.raster) <- paste(RASTER.species.name, "_", stack.title , "_zones", sep="")
    raster::writeRaster(working.raster, filename=rasterfull, progress='text', overwrite=TRUE, format=RASTER.format, datatype=RASTER.datatype, NAflag=RASTER.NAflag)
#
    if (KML.out == T) {
        nc <- nrow(centroid.object$centroids)
        raster::KML(working.raster, filename=kmlfull, col = grDevices::rainbow(n = nc, start = 0.2, end = 0.8), colNA = 0, 
            blur=KML.blur, maxpixels=KML.maxpixels, overwrite=TRUE, breaks = c(0:nc))
    }

    cat(paste("\n", "zones provided in folder: ", getwd(), "//ensembles//zones", "\n", sep=""))
    zones.raster <- raster::raster(rasterfull)
    return(zones.raster)
}

