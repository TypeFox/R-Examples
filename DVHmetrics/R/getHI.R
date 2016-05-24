## we need ... because getMetric() will also pass parameters
## intended for other functions through ...
getHI <-
function(x, ...) {
    UseMethod("getHI")
}

getHI.DVHs <-
function(x, ...) {
    Dmetr <- getMetric(x, metric=c("D2%", "D50%", "D98%"))
    HI <- (Dmetr[["D2%"]] - Dmetr[["D98%"]]) / Dmetr[["D50%"]]

    data.frame(HI=HI,
               patID=x$patID,
               structure=x$structure,
               stringsAsFactors=FALSE)
}

getHI.DVHLst <-
function(x, EUDa, EUDfd=NULL, EUDab=NULL, ...) {
    HIl  <- Map(getHI, x)
    HIdf <- do.call("rbind", HIl)
    rownames(HIdf) <- NULL
    HIdf
}

getHI.DVHLstLst <-
function(x, EUDa, EUDfd=NULL, EUDab=NULL, ...) {
    HIl  <- Map(getHI, x)
    HIdf <- do.call("rbind", HIl)
    rownames(HIdf) <- NULL
    HIdf
}
