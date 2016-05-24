`mefaCrosstab` <-
function(x, segment=FALSE, nested=FALSE, drop.zero=FALSE)
{
    if (!is.stcs(x))
        stop("'x' must be of class 'stcs'")
#    ss <- stcs(x, expand = !attr(x, "expand"), drop.zero=drop.zero)
    ss <- stcs(x, expand = FALSE, drop.zero=drop.zero)
    if (!segment) {
 #       out <- as.matrix(table(ss$samp, ss$taxa))
        out <- xtabs(count ~ samp + taxa, ss)
        class(out) <- NULL
        if (attr(ss, "zero.count") && !drop.zero)
            out <- out[, -which(colnames(out) %in% attr(ss, "zero.pseudo"))]}
    if (segment) {
        out <- list()
        nsegm <- if (attr(ss, "zero.count"))
            nlevels(ss$segm) - 1 else nlevels(ss$segm)
        if (nsegm == 1) {
            out <- NULL
            } else {
            for (i in 1:nsegm) {
                if (drop.zero || (levels(ss$segm)[i] != attr(ss, "zero.pseudo"))) {
#                    out[[i]] <- as.matrix(table(ss$samp, ss$tax, ss$segm)[,,i])
                    out[[i]] <- xtabs(count ~ samp + taxa + segm, ss)[,,i]
                    class(out[[i]]) <- NULL
                if (attr(ss, "zero.count") && !drop.zero)
                    out[[i]] <- out[[i]][, -which(colnames(out[[i]]) %in% attr(ss, "zero.pseudo"))]
                names(out)[i] <- levels(ss$segm)[i]}}
            if (nested) {
                lnam <- names(out)
                for (i in 2:length(out)) {
                    out[[i]] <- out[[(i-1)]] + out[[i]]
                    names(out)[[i]] <- paste(lnam[1], "-", lnam[i], sep="")}}
            }
        }
return(out)
}

