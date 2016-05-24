##############################################################################
## package 'secr'
## addCovariates.R
## 2011-11-01
## 2013-01-19 handles missing values in character fields
## 2014-08-05 strict argument
## 2014-08-05 relax requirement for object to be traps or mask:
## 2014-08-05 may now be any numeric vector that can be formed into a 2-column matrix
###############################################################################

addCovariates <- function (object, spatialdata, columns = NULL, strict = FALSE) {
    if (!(inherits(object, 'mask') | inherits(object, 'traps')))
        ## stop ("require mask or traps object")
        object <- matrix(unlist(object), ncol = 2)
    if (!ms(object) & ms(spatialdata))
        stop ("mismatch of single session object, multisession spatialdata")

    if (ms(object)) {
        ## allow multiple sessions, and session-specific data sources
        nsession <- length(object)
        out <- object
        for (session in 1:nsession) {
            if (ms(spatialdata))
                out[[session]] <- addCovariates(out[[session]], spatialdata[[session]])
            else
                out[[session]] <- addCovariates(out[[session]], spatialdata)
        }
        out
    }
    else {

        if (is.character(spatialdata))
            type <- "shapefile"
        else if (inherits(spatialdata, "SpatialPolygonsDataFrame"))
            type <- "SPDF"
        else if (inherits(spatialdata, "SpatialGridDataFrame"))
            type <- "SGDF"
        else if (inherits(spatialdata, "mask"))
            type <- "mask"
        else if (inherits(spatialdata, "traps"))
            type <- "traps"
        else
            stop ("spatialdata type unrecognised or unsupported")

        if (type == "shapefile") {
            polyfilename <- spatialdata  ## strip shp?
            if (!requireNamespace('maptools', quietly = TRUE))
                stop("requires maptools")
            spatialdata <- maptools::readShapePoly(polyfilename)
        }
        if (type %in% c("shapefile", "SPDF", "SGDF")) {
            xy <- matrix(unlist(object), ncol = 2)
            xy <- SpatialPoints(xy)
            ## 2014-12-11
            proj4string(spatialdata) <- CRS()
            df <- sp::over (xy, spatialdata)
        }
        else {
            ## nearest point algorithm
            if (is.null(covariates(spatialdata)))
                stop ("spatialdata does not have covariates")
            index <- nearesttrap(object, spatialdata)
            df <- covariates(spatialdata)[index,, drop=FALSE]
            ## new argument 2014-08-05
            if (strict) {
                incell <- function (xy, m, mask) {
                    sp2 <- spacing(mask) / 2
                    mxy <- mask[m,]
                    ((xy[,1] + sp2) >= mxy[,1]) &
                    ((xy[,1] - sp2) <= mxy[,1]) &
                    ((xy[,2] + sp2) >= mxy[,2]) &
                    ((xy[,2] - sp2) <= mxy[,2])

                }
                cellOK <- incell(object, index, spatialdata)
                df[!cellOK,] <- NA
                if (any(!cellOK))
                    warning ("some requested points lie outside mask")
            }
        }

        ## select requested columns
        if (!is.null(columns))
            df <- df[,columns, drop = FALSE]

        ## check new covariates OK
        fn <- function(x) {
            if (is.numeric(x))
                !any(is.na(x))
            else
                !any((x == "") | is.na(x))
        }
        OK <- all(apply(df, 2, fn))
        if (!OK)
            warning ("missing values among new covariates")

        ## insert new covariates and return object
        rownames(df) <- 1:nrow(df)
        if (is.null(covariates(object)))
            covariates(object) <- df
        else
            covariates(object) <- cbind(covariates(object), df)
        object
    }
}

