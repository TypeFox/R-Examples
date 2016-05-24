
".derivStats" <- function(x, diveNo)
{
    ## Value: Matrix with one row per dive, keeping order in 'x'
    ## --------------------------------------------------------------------
    ## Arguments: x=TDRcalibrate object; diveNo=numeric vector specifying
    ## which dives to obtain derivative statistics.
    ## --------------------------------------------------------------------
    ## Purpose: Provide summary statistics (summary() and sd()) of
    ## derivatives for descent, bottom, and ascent phases for each dive.
    ## --------------------------------------------------------------------
    ## Author: Sebastian P. Luque
    ## --------------------------------------------------------------------
    if (missing(diveNo)) {
        diveNo <- seq(max(getDAct(x, "dive.id")))
    } else {
        diveNo <- sort(unique(diveNo))
    }
    summarize.phase <- function(diveNo, phase, label) {
        der <- getDiveDeriv(x, diveNo=diveNo, phase=phase)
        der.summ <- summary(der$y)
        names(der.summ) <- gsub("[ \\.]", "", tolower(names(der.summ)))
        der.sd <- sd(der$y)
        der.m <- c(der.summ, sd=der.sd)
        names(der.m) <- paste(label, names(der.m), sep=".")
        der.m
    }
    d <- do.call(rbind,
                 lapply(diveNo, summarize.phase, "descent", "descD"))
    b <- do.call(rbind,
                 lapply(diveNo, summarize.phase, "bottom", "bottD"))
    a <- do.call(rbind,
                 lapply(diveNo, summarize.phase, "ascent", "ascD"))
    cbind(d, b, a)
}

"diveStats" <- function(x, depth.deriv=TRUE)
{
    ## Value: A data frame with per-dive statistics
    ## --------------------------------------------------------------------
    ## Arguments: x=object of class TDRcalibrate
    ## --------------------------------------------------------------------
    ## Author: Sebastian Luque
    ## --------------------------------------------------------------------
    if (!is(x, "TDRcalibrate")) stop("x must be a TDRcalibrate object")
    zvtdr <- getTDR(x)                     # fully calibrated object
    interval <- getDtime(zvtdr)            # sampling interval
    diveid <- getDAct(x, "dive.id")        # dive IDs
    postdiveid <- getDAct(x, "postdive.id")          # postdive IDs
    ok <- which(diveid > 0 & diveid %in% postdiveid) # diving subscripts
    dtimes <- getTime(zvtdr)[ok]                     # diving times
    ddepths <- getDepth(zvtdr)[ok]                   # diving depths
    dids <- diveid[ok]                               # dive IDs
    dphases <- getDPhaseLab(x)[ok]                   # dive phase labels
    okpd <- which(postdiveid %in% unique(dids)) # postdive subscripts
    pdtimes <- getTime(zvtdr)[okpd]          # required postdive times
    pddepths <- getDepth(zvtdr)[okpd]        # required postdive depths
    pdids <- postdiveid[okpd]                # required postdive IDs

    postdive.dur <- tapply(pdtimes, pdids, function(k) {
        difftime(k[length(k)], k[1], units="secs")
    })

    dtimestz <- attr(dtimes, "tzone")
    if (!is(zvtdr, "TDRspeed")) {
        td <- data.frame(dphases, dtimes, ddepths)
        perdive <- do.call(rbind, by(td, dids, oneDiveStats, interval=interval))
        res <- data.frame(perdive, postdive.dur)
        for (i in 1:3) res[, i] <- .POSIXct(res[, i], dtimestz)
    } else {
        dspeeds <- getSpeed(zvtdr)[ok]  # diving speeds
        td <- data.frame(dphases, dtimes, ddepths, dspeeds)
        perdive <- do.call(rbind, by(td, dids, oneDiveStats, interval=interval,
                                     speed=TRUE))
        ## for postdive total distance and mean speed
        ptd <- matrix(c(pdtimes, getSpeed(zvtdr)[okpd]), ncol=2)
        pdv <- do.call(rbind, by(ptd, pdids, .speedStats))
        res <- data.frame(perdive, postdive.dur, postdive.tdist=pdv[, 1],
                          postdive.mean.speed=pdv[, 2], row.names=NULL)
        for (i in 1:3) res[, i] <- .POSIXct(res[, i], dtimestz)
    }

    if (depth.deriv) {
        data.frame(res, .derivStats(x, diveNo=dids))
    } else res
}
