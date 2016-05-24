RCS <-
function(data.tc, level=.01) {
# radiometric control sets, HSNG199
# takes output of tasscap
#    list(Brightness=Brightness, Greenness=Greenness, Wetness=Wetness)
# or named data frame, which is also a list

    rcsgrid <- data.tc$Brightness
    
    brightness <- as.vector(as.matrix(data.tc$Brightness))
    greenness <- as.vector(as.matrix(data.tc$Greenness))
        
    bright.llevel <- quantile(brightness, level, na.rm=TRUE)
    bright.ulevel <- quantile(brightness, 1-level, na.rm=TRUE)
    green.level <- quantile(greenness, level, na.rm=TRUE)
    
    rcsmask <- ifelse(brightness < bright.llevel & greenness < green.level, 1, 0)
    rcsmask <- ifelse(brightness > bright.ulevel & greenness < green.level, 1, rcsmask)

    # return the same structure as the input values
    if(class(rcsgrid) == "SpatialGridDataFrame") 
        rcsgrid@data[,1] <- rcsmask
    else if(is.data.frame(rcsgrid)) 
        rcsgrid <- data.frame(matrix(rcsmask, nrow=nrow(rcsgrid), ncol=ncol(rcsgrid)))
    else if(is.matrix(rcsgrid)) 
        rcsgrid <- matrix(rcsmask, nrow=nrow(rcsgrid), ncol=ncol(rcsgrid))
    else # return a vector 
        rcsgrid <- rcsmask

    rcsgrid
}

