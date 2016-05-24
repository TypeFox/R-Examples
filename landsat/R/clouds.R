
clouds <- function(band1, band6, level = 0.0014, buffer=5) {
    # simple function to create a cloud mask
    # uses bands 6 and 1 as specified in the arguments
    # cloud target: band6 is low; band1 is high
   
    # clouds have high reflectance in band1 and
    # low temperature (band 6)
    # thus ratio band1/band6 is high for clouds

    if(is.character(band1)) {
        band1 <- read.asciigrid(band1)
        results <- band1
        dims <- band1@grid@cells.dim
        band1 <- band1@data[,1]
    } else {
        if(class(band1) == "SpatialGridDataFrame") {
            results <- band1
            dims <- band1@grid@cells.dim
            band1 <- band1@data[,1]
        } else {
            results <- band1
            dims <- dim(as.matrix(band1))
            band1 <- as.vector(as.matrix(band1))
        }
    } 
    
    if(is.character(band6)) {
        band6 <- read.asciigrid(band6)
        band6 <- band6@data[,1]
    } else {
        if(class(band6) == "SpatialGridDataFrame") {
            band6 <- band6@data[,1]
        } else {
            band6 <- as.vector(as.matrix(band6))
        }
    } 

           
    cloudmask <- ifelse(band1/band6 > level, 1, 0) # 6 low AND 1 high
    ## want to add a buffer around identified areas
    cloudmask <- movingwindow(matrix(cloudmask,  nrow=dims[1], ncol=dims[2]), matrix(rep(1, buffer*buffer), buffer, buffer))
    cloudmask <- ifelse(as.vector(cloudmask) > 0, 1, NA)

    # return the same structure as the input values
    if(class(results) == "SpatialGridDataFrame")
        results@data[,1] <- cloudmask
    else if(is.data.frame(results))
        results <- data.frame(matrix(cloudmask, nrow=nrow(results), ncol=ncol(results)))
    else if(is.matrix(results))
        results <- matrix(cloudmask, nrow=nrow(results), ncol=ncol(results))
    else # return a vector 
        results <- cloudmask
    
    results
}

