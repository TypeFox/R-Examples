mci <- function (x) 
{
    syncDat <- split(x, x$syncID)
    if (length(unique(sapply(syncDat, nrow))) > 1) {
        stop("not all synchronization events have a record for each
             individual")
    }
    syncID <- unique(x$syncID)

    ## determine pairs of subsequent synchronized events and
    ## get lower indexes of shift
    shift.lo <- which(c(diff(syncID) == 1, FALSE))

    ## if no subsequent synchronized events present end function
    if (!length(shift.lo)) {
        return(message("no subsequent synchronization events found"))
    }

    ## get higher indexes of shift
    shift.hi <- shift.lo + 1

    ## calculate shift in the x coordinate
    x.shift <- mapply(function(f1, f2) syncDat[[f2]]$utm.easting - 
        syncDat[[f1]]$utm.easting, f1 = shift.lo, f2 = shift.hi)

    ## calculate shift in the y coordinate
    y.shift <- mapply(function(f1, f2) syncDat[[f2]]$utm.northing - 
        syncDat[[f1]]$utm.northing, f1 = shift.lo, f2 = shift.hi)

    ## calculate MCI indexes
    mci <- 1 - 0.5 * (apply(x.shift, 2, function(f) sum(abs(f - 
        mean(f)))/sum(abs(f))) + apply(y.shift, 2, function(f) sum(abs(f - 
        mean(f)))/sum(abs(f))))
    return(mci)
}
