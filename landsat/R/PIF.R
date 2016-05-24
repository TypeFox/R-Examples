PIF <-
function(band3, band4, band7, level=.99) {
# identify pseudo-invariant features after SSV1988
   
    if(is.character(band3)) {
        band3 <- read.asciigrid(band3)
        pifgrid <- band3
        band3 <- band3@data[,1]
    } else {
        pifgrid <- band3
        band3 <- as.vector(as.matrix(band3))
    } 
    
    if(is.character(band4)) {
        band4 <- read.asciigrid(band4)@data[,1]
    } else {
        band4 <- as.vector(as.matrix(band4))
    }

    if(is.character(band7)) {
        band7 <- read.asciigrid(band7)@data[,1]
    } else {
        band7 <- as.vector(as.matrix(band7))
    }

    band43 <- band4/band3
        
    band43.level <- quantile(band43, 1-level, na.rm=TRUE)
    band7.level <- quantile(band7, level, na.rm=TRUE)
    
    pifmask <- ifelse(band43 < band43.level & band7 > band7.level & band7 < 255, 1, 0)
    
    # return the same structure as the input values
    if(class(pifgrid) == "SpatialGridDataFrame")
        pifgrid@data[,1] <- pifmask
    else if(is.data.frame(pifgrid))
        pifgrid <- data.frame(matrix(pifmask, nrow=nrow(pifgrid), ncol=ncol(pifgrid)))
    else if(is.matrix(pifgrid))
        pifgrid <- matrix(pifmask, nrow=nrow(pifgrid), ncol=ncol(pifgrid))
    else # return a vector 
        pifgrid <- pifmask
    
    pifgrid
}

