
".depthFilter" <- function(depth, k, probs, depth.bounds, na.rm)
{
    ## Value: A matrix with the filtered depths at each step, and corrected
    ## depth as last column, so it has dimensions length(depth) rows by
    ## length(k) + 1 columns
    ## --------------------------------------------------------------------
    ## Arguments: depth=depth vector; k=vector of window width integers to
    ## be applied sequentially (coerced to integer in runquantile());
    ## probs=vector of quantiles to extract at each step indicated by k (so
    ## must be as long as k); depth.bounds=minimum and maximum depth bounds
    ## where search should be restricted to; na.rm=do we remove NA from
    ## depth before filtering (recommended if there are no level shifts)?
    ## --------------------------------------------------------------------
    ## Purpose: Calculate running quantiles on sequential filters, starting
    ## with the original depth vector and correct depth as: original depth
    ## - the filtered/smoothed depth
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    if (length(k) != length(probs))
        stop("k and probs should have the same length")
    if (length(depth.bounds) != 2 && !is.numeric(depth.bounds))
        stop("depth.bounds must be a 2-element numeric vector")
    d.na <- is.na(depth)
    d.in.bounds <- depth > depth.bounds[1] & depth < depth.bounds[2]
    if (na.rm) {
        d.ok <- !d.na & d.in.bounds     # logical
    } else d.ok <- which(d.in.bounds | is.na(depth))   # numeric
    filters <- matrix(depth, ncol=1)
    for (i in seq(length(k))) {
        filters <- cbind(filters, depth)
        dd <- filters[d.ok, i]
        filters[d.ok, i + 1] <- caTools::runquantile(dd, k=k[i],
                                                     probs=probs[i])
        ## Linear interpolation for depths out of bounds
        offbounds <- which(!d.in.bounds)
        offbounds.fun <- approxfun(seq(length(depth))[d.in.bounds],
                                   filters[d.in.bounds, i + 1], rule=2)
        filters[offbounds, i + 1] <- offbounds.fun(offbounds)                                   
        ## NA input should be NA output regardless of na.rm
        filters[d.na, i + 1] <- NA
    }
    dnames <- c("depth",
                paste("smooth", paste(k, probs, sep="_"), sep="."))
    dimnames(filters)[[2]] <- dnames
    depth.adj <- depth - filters[, ncol(filters)]
    cbind(filters[, -1], depth.adj)     # don't output original depth
}


".zoc" <- function(time, depth, method, control)
{
    ## Value: Vector with all corrected depths for time windows based on
    ## 'method'. Make depths < 0 equal 0
    ## --------------------------------------------------------------------
    ## Arguments: time=POSIXct; depth=uncalibrated depth; method=method to
    ## use for adjusting depth; control=named list of control parameters to
    ## be used for the chosen method (zoc.filter: k, probs, na.rm; offset:
    ## offset)
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque (thanks to Roland Fried for suggesting
    ## something like method="filter")
    ## --------------------------------------------------------------------
    switch(method,
           visual = {zoclims <- plotTDR(time, depth)
                     dev.off()
                     if (length(zoclims) == 0) {
                         message("No ZOC performed.")
                     } else {
                         ## Subtract the ZOC from depths in each window
                         ## If there's any overlap, the latest window prevails
                         zocdives <- lapply(zoclims, function(x) {
                             newdep <- depth[time >= x[[1]][1] &
                                             time < x[[1]][2]] - x[[2]][1]
                             indsnewdep <- which(time >= x[[1]][1] &
                                                 time < x[[1]][2])
                             cbind(indsnewdep, newdep)
                         })
                         message(paste("ZOC procedure completed\n",
                                       length(zocdives),
                                       " time windows have been corrected",
                                       sep=""))
                         zocdives <- do.call("rbind", zocdives)
                         ## Correct depths within chosen time windows
                         depth[zocdives[, 1]] <- zocdives[, 2]
                     }},
           offset = {offset <- control$offset
                     depth <- depth - offset},
           filter = {k <- control$k
                     probs <- control$probs
                     depth.bounds <- control$depth.bounds
                     na.rm <- control$na.rm
                     depthmtx <- .depthFilter(depth, k, probs,
                                              depth.bounds, na.rm)
                     depth <- depthmtx[, ncol(depthmtx)]})
    ## Turn all negative and NA depths into zeroes (we don't care now about
    ## dry time, since we've already done that, and we need surface and dry
    ## time to be zero from this point on).
    depth[depth < 0 & !is.na(depth)] <- 0
    depth
}
