`checkZone` <-
function(from, zones, checkz, nofsets, ib, ellipse) {
    
    fromz <- ib[ib$s == nofsets & ib$v == as.numeric(ellipse) & ib$i == from, ]
    toz <- ib[ib$s == nofsets & ib$v == as.numeric(ellipse) & ib$i %in% zones[!checkz], ]
    toz <- toz[toz$b %in% fromz$b, , drop = FALSE]
    
    if (nrow(toz) > 0) {
        zs <- sort(unique(toz$i))
        
        checkz[as.character(zs)] <- TRUE
        for (i in zs) {
            checkz <- checkz | Recall(i, zones, checkz, nofsets, ib, ellipse)
        }
    }
    
    return(checkz)
}

