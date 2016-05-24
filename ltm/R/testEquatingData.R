testEquatingData <-
function (DataList, AnchoringItems = NULL) {
    if (!is.list(DataList) || any(sapply(DataList, function (x) !inherits(x, "matrix") && !inherits(x, "data.frame"))))
        stop("\n'DataList' must be a list containing either data.frames or matrices.")
    dat.lis <- lapply(DataList, data.matrix)
    if (!is.null(AnchoringItems)) {
        if (!inherits(AnchoringItems, "matrix") && !inherits(AnchoringItems, "data.frame"))
            stop("\n'AnchoringItems' must be either a data.frame or a matrix.")
        AnchoringItems <- data.matrix(AnchoringItems)
        if (is.null(colnames(AnchoringItems)))
            colnames(AnchoringItems) <- paste("AnchrItm", 1:ncol(AnchoringItems), sep = "")
        if (any(ind <- colnames(AnchoringItems) %in% unlist(lapply(dat.lis, colnames))))
            colnames(AnchoringItems)[ind] <- paste("AnchrItm", colnames(AnchoringItems)[ind], sep = ".")
        dat.lis <- lapply(dat.lis, function (x) cbind(AnchoringItems, x))
    }
    itms.nams <- unlist(lapply(dat.lis, colnames))
    uniq.itms <- unique(itms.nams)
    nrs <- sapply(dat.lis, NROW)
    n.max <- sum(nrs)
    p.max <- length(uniq.itms)
    out <- matrix(NA, n.max, p.max)
    colnames(out) <- uniq.itms
    rownames(out) <- unlist(lapply(dat.lis, rownames))
    ind1 <- c(1, cumsum(nrs[-length(nrs)]) + 1)
    ind2 <- cumsum(nrs)
    for (i in seq(dat.lis)) {
        dat <- dat.lis[[i]]
        nam.ind <- uniq.itms %in% colnames(dat)
        out[seq(ind1[i], ind2[i]), nam.ind] <- dat[, uniq.itms[nam.ind]]
    }
    out <- out[, order(colSums(is.na(out)))]
    ind <- sum(colSums(is.na(out)) == 0)
    attr(out, "items") <- if (!is.null(AnchoringItems)) {
        rep(c(" *", ""), c(ncol(AnchoringItems), ncol(out) - ncol(AnchoringItems)))
    } else {
        rep(c(" *", ""), c(ind, p.max - ind))
    }
    attr(out, "anchoring") <- !is.null(AnchoringItems)
    out
}
