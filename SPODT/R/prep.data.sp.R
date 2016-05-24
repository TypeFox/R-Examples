prep.data.sp <-
function(data)
{
    z <- as.vector(by(data$z, data$loc, sum))
    n <- as.vector(table(data$loc))
    
    data.sp <- unique(cbind.data.frame(data$loc, data$x, data$y))    
    colnames(data.sp) <- c("loc","x","y")    
    data.sp <- rbind(data.sp)[order(data.sp$loc),]
    data.sp <- rbind(cbind.data.frame(data.sp$loc, z, n))[order(data.sp$x, data.sp$y),]
        
    colnames(data.sp) <- c("loc","z","n")
    rownames(data.sp) <- rownames(data.sp[which(data.sp$loc %in% data$loc),])

    return(data.sp)
}
