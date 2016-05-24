
##interGC <- function(x, dur = NULL, intfun = NULL) {
##    if (is.null(intfun)) {
##        require(geosphere)
##        intfun <- gcIntermediate
##    }

##    tor <- getTORnames(x)
##    ids <- x[[tor[2]]]
##    times <- x[[tor[1]]]
##    if (is.null(dur))  stop("dur must be specified for minimum interpolation duration")
    ##dur <- min(sapply(split(times, ids), function(x) min(diff(unclass(x)))))

##    triplist <- split(x, ids)
##    newtriplist <- vector("list", length(triplist))
##    for (itrip in seq_along(triplist)) {
##        icoords <- coordinates(triplist[[itrip]])
 ##       itimes <- triplist[[itrip]][[tor[1]]]
 ##       newcoords <- vector("list", nrow(icoords) - 1)
 ##       dntime <- pmax(ceiling(diff(unclass(itimes)) / dur), 3)
 ##       for (ipt in seq_len(nrow(icoords)-1)) {
 ##           newcoords[[ipt]] <- data.frame(intfun(icoords[ipt, ], icoords[ipt+1,], n = dntime[ipt] - 2, addStartEnd = TRUE),
 ##                              gmt = seq(itimes[ipt], itimes[ipt+1], length = dntime[ipt]), id = rep(ids[itrip], dntime[ipt]))

   ##     }
    ##    newtriplist[[itrip]] <- do.call(rbind, newcoords)
   ## }
   ## newtrip <- do.call(rbind, newtriplist)
   ## names(newtrip) <- c(colnames(icoords), tor)
   ## newtrip <- newtrip[!duplicated(newtrip), ]
   ## coordinates(newtrip) <- 1:2
   ## proj4string(newtrip) <- CRS(proj4string(x))
   ## trip(newtrip, tor)
##}


