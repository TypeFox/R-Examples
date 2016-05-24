imageplot.bma <-
function (bma.out, color = c("red", "blue", "#FFFFD5"), order = c("input", 
    "probne0", "mds"), ...) 
{
    clr <- color
    if (length(color) == 1) {
        if (color == "default") 
            clr <- c("#FF0000", "#FF0000", "#FFFFD5")
        if (color == "blackandwhite") 
            clr <- c("black", "black", "white")
    }
    keep.mar <- par(mar = c(5, 6, 4, 2) + 0.1)
    nmodel <- nrow(bma.out$which)
    which <- bma.out$which
    probne0 <- bma.out$probne0
    if (class(bma.out) == "bic.surv") 
        mle <- bma.out$mle
    else mle <- bma.out$mle[, -1, drop = FALSE]
    nvar <- ncol(mle)
    rownms <- bma.out$namesx
    if (ifelse(!is.null(bma.out$factor.type), bma.out$factor.type, 
        FALSE)) {
        which <- matrix(NA, ncol = nvar, nrow = nmodel)
        probne0 <- rep(NA, times = nvar)
        rownms <- rep(NA, times = nvar)
        assign <- bma.out$assign
        offset <- 1
        if (class(bma.out) == "bic.surv") 
            offset <- 0
        assign[[1]] <- NULL
        for (i in 1:length(assign)) {
            probne0[assign[[i]] - offset] <- bma.out$probne0[i]
            which[, assign[[i]] - offset] <- bma.out$which[, 
                i]
            nm <- names(bma.out$output.names)[i]
            if (!is.na(bma.out$output.names[[i]][1])) 
                nm <- paste(nm, bma.out$output.names[[i]][-1], 
                  sep = ".")
            rownms[assign[[i]] - offset] <- nm
        }
    }
    ordr.type <- match.arg(order)
    if (ordr.type == "probne0") 
        ordr <- order(-probne0)
    else if (ordr.type == "mds") {
        postprob.rep <- matrix(bma.out$postprob, ncol = nvar, 
            nrow = nmodel)
        k11 <- t(which + 0) %*% ((which + 0) * postprob.rep)
        k00 <- t(1 - which) %*% ((1 - which) * postprob.rep)
        k01 <- t(which + 0) %*% ((1 - which) * postprob.rep)
        k10 <- t(1 - which) %*% ((0 + which) * postprob.rep)
        ktau <- 4 * (k00 * k11 - k01 * k10)
        dissm <- 1 - abs(ktau)
        diag(dissm) <- 0
        ordr <- order(as.vector(cmdscale(dissm, k = 1)))
    }
    else ordr <- 1:nvar
    ordr <- rev(ordr)
    postprob <- bma.out$postprob
    which <- which[, ordr, drop = FALSE]
    mle <- mle[, ordr, drop = FALSE]
    rownms <- rownms[ordr]
    color.matrix <- (which) * (2 - (mle > 0)) + 3 * (!which)
    par(las = 1)
    image(c(0, cumsum(postprob)), 1:nvar, color.matrix, col = clr, 
        xlab = "Model #", ylab = "", xaxt = "n", yaxt = "n", 
        xlim = c(0, 1), main = "Models selected by BMA", ...)
    xat <- (cumsum(postprob) + c(0, cumsum(postprob[-nmodel])))/2
    axis(1, at = xat, labels = 1:nmodel, ...)
    axis(2, at = 1:nvar, labels = rownms, ...)
    par(mar = keep.mar)
}

