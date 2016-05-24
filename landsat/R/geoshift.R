geoshift <-
function(mat, padx, pady, shiftx, shifty, nodata=NA)
{

    results <- mat

    if(is.data.frame(mat)) mat <- as.matrix(mat)
    if(is.matrix(mat)) {
        # pad a matrix and slide it over a given amount
        newmat <- matrix(nodata, nrow=(nrow(mat) + 2 * padx), ncol=(ncol(mat) + 2 * pady))
        newmat[(shiftx + padx + 1):(nrow(mat) + shiftx + padx), (shifty + pady + 1):(ncol(mat) + shifty + pady)] <- mat
    }
    if(class(mat) == "SpatialGridDataFrame") {
        mat.data <- as.matrix(mat)
        mat.data <- geoshift(mat.data, padx, pady, shiftx, shifty, nodata=nodata)
        mat@data <- data.frame(as.vector(mat.data))

        # S4 class
        mat.grid <- mat@grid
        mat.grid@cellcentre.offset[1] <- mat.grid@cellcentre.offset[1] - (padx * mat.grid@cellsize[1])
        mat.grid@cellcentre.offset[2] <- mat.grid@cellcentre.offset[2] - (pady * mat.grid@cellsize[2])
        mat.grid@cells.dim[1] <- as.integer(mat.grid@cells.dim[1] + 2*padx)
        mat.grid@cells.dim[2] <- as.integer(mat.grid@cells.dim[2] + 2*pady)
        mat@grid <- mat.grid

## obsoleted by sp_0.9-94
#-#        # not S4 class
#-#        mat.coords <- coordinates(mat)
#-#        mat.coords[1, 1] <- mat.coords[1, 1] - (padx * mat.grid@cellsize[1])
#-#        mat.coords[2, 1] <- mat.coords[2, 1] + (padx * mat.grid@cellsize[1])
#-#        mat.coords[1, 2] <- mat.coords[1, 2] - (pady * mat.grid@cellsize[2])
#-#        mat.coords[2, 2] <- mat.coords[2, 2] + (pady * mat.grid@cellsize[2])
#-#        mat@coords <- mat.coords
         
        # not S4 class
        mat.bbox <- bbox(mat)
        mat.bbox[1, "min"] <- mat.bbox[1, "min"] - (padx * mat.grid@cellsize[1])
        mat.bbox[1, "max"] <- mat.bbox[1, "max"] + (padx * mat.grid@cellsize[1])
        mat.bbox[2, "min"] <- mat.bbox[2, "min"] - (pady * mat.grid@cellsize[2])
        mat.bbox[2, "max"] <- mat.bbox[2, "max"] + (pady * mat.grid@cellsize[2])
        mat@bbox <- mat.bbox
    }
    
    # return the same structure as the input values
    if(class(results) == "SpatialGridDataFrame")
        results <- mat
    else if(is.data.frame(results))
        results <- data.frame(newmat)
    else if(is.matrix(results))
        results <- newmat
    else # no result
        results <- NA
    
    results

}

