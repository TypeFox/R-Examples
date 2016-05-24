# get the count number of time pairs output by GetRawCov
# Output: a data.frame of three columns: t1, t2, count
GetCount <- function(tpairs) {
    # browser()
    tab <- table(tpairs[, 1], tpairs[, 2])
    pts <- sort(unique(as.numeric(tpairs)))
    ret <- data.frame(expand.grid(pts, pts), as.numeric(tab))
    names(ret) <- c('t1', 't2', 'count')
    ret <- ret[ret$count != 0, ]
    
    return(ret)
}



