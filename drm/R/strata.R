"strata" <-
function (..., na.group = FALSE, shortlabel = FALSE) 
{
    words <- as.character((match.call())[-1])
    if (!missing(na.group)) 
        words <- words[-1]
    allf <- list(...)
    if (length(allf) == 1 && is.list(ttt <- unclass(allf[[1]]))) {
        allf <- ttt
        words <- names(ttt)
    }
    nterms <- length(allf)
    what <- allf[[1]]
    if (is.null(levels(what))) 
        what <- factor(what)
    levs <- unclass(what) - 1
    wlab <- levels(what)
    if (na.group && any(is.na(what))) {
        levs[is.na(levs)] <- length(wlab)
        wlab <- c(wlab, "NA")
    }
    if (shortlabel) 
        labs <- wlab
    else labs <- paste(words[1], wlab, sep = "=")
    for (i in (1:nterms)[-1]) {
        what <- allf[[i]]
        if (is.null(levels(what))) 
            what <- factor(what)
        wlab <- levels(what)
        wlev <- unclass(what) - 1
        if (na.group && any(is.na(wlev))) {
            wlev[is.na(wlev)] <- length(wlab)
            wlab <- c(wlab, "NA")
        }
        if (!shortlabel) 
            wlab <- format(paste(words[i], wlab, sep = "="))
        levs <- wlev + levs * (length(wlab))
        labs <- paste(rep(labs, rep(length(wlab), length(labs))), 
            rep(wlab, length(labs)), sep = ", ")
    }
    levs <- levs + 1
    ulevs <- sort(unique(levs[!is.na(levs)]))
    levs <- match(levs, ulevs)
    labs <- labs[ulevs]
    factor(levs, labels = labs)
}
