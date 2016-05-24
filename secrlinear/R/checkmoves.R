############################################################################################
## package 'secrlinear'
## checkmoves.R
## last changed 2014-12-01
############################################################################################

checkmoves <- function (CH, accept = c(0, 1000), userdist, mask, showall = TRUE, silent = FALSE) {
    if (ms(CH))
        stop("not ready for multi-session data")
    trps <- traps(CH)
    if (missing(userdist)) {
        if (missing(mask))
            userdist <- as.matrix(dist(trps))
        else {
            if (!inherits(mask, 'linearmask'))
                stop ("linearmask expected")
            userdist <- networkdistance(trps, trps, mask)
        }
    }
    if (any(dim(userdist) != nrow(trps)))
        stop("userdist should be matrix of distances between detectors")
    rownames(userdist) <- rownames(trps)
    m <- moves(CH, userdist = userdist)
    ## infmoves <- sapply(m, function(x) any(!is.finite(x)))
    badmoves <- sapply(m, function(x) any((x < accept[1]) | (x > accept[2])))
    out <- list(badmoves = badmoves)
    if (!any(badmoves)) {
        if (!silent) 
            cat('All OK\n')
    }
    else {
        if (!silent) 
            cat(sum(badmoves), ' animal(s) with movement(s) outside range ', 
                accept[1], ' to ', accept[2], '\n')
        out$CH <- suppressWarnings(subset(CH, names(m)[badmoves]))
        trap <- trap(out$CH)
        move <- m[badmoves]
        move <- round(unlist(lapply(move, c, NA)), 2)
        
        df <- data.frame(
            ID = animalID(out$CH),
            occasion = occasion(out$CH),
            detector = trap,
            stringsAsFactors = FALSE
            )
        if (!missing(mask))
            df$LineID <- covariates(mask)$LineID[nearesttrap(subset(traps(CH),trap), mask)]        
        out$df <- df[order(df$ID, df$occasion, df$detector),]
        rownum <- match(out$df$detector, rownames(userdist))
        out$df$movement <- c(userdist[cbind(rownum[-nrow(out$df)], rownum[-1])], NA)
        out$df$movement[ match(unique(out$df$ID), out$df$ID) - 1] <- NA
        out$df$movement <- round(out$df$movement, 2)
        if (!showall) {
            OK <- ((out$df$movement < accept[1]) | (out$df$movement >  accept[2])) & !is.na(out$df$movement)
            out$df <- out$df[OK,]
        }
        rownames(out$df) <- 1:nrow(out$df)
    }
    invisible(out)
}

