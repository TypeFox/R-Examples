OBsProb <-
function (X, y, abeta=1, bbeta=1, blk , mFac, mInt, nTop) 
{
    X <- as.matrix(X)
    y <- unlist(y)
    if (length(y) != nrow(X)) 
        stop("X and y should have the same number of observations")
    if (blk == 0) {
        ifelse(is.null(colnames(X)), faclab <- paste("F", seq(ncol(X)), 
            sep = ""), faclab <- colnames(X))
        colnames(X) <- faclab
    }
    else {
        if (is.null(colnames(X))) {
            faclab <- paste("F", seq(ncol(X) - blk), sep = "")
            blklab <- paste("B", seq(blk), sep = "")
            colnames(X) <- c(blklab, faclab)
        }
        else {
            faclab <- colnames(X)[-seq(blk)]
            blklab <- colnames(X)[seq(blk)]
        }
    }
    rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "r")
    storage.mode(X) <- "double"
    Y <- as.double(y)
    N <- as.integer(nrow(X))
    COLS <- as.integer(ncol(X) - blk)
    BLKS <- as.integer(blk)
    MXFAC <- as.integer(max(1, mFac))
    MXINT <- as.integer(mInt)
    NTOP <- as.integer(nTop)
    mdcnt <- as.integer(0)
    abeta <- as.integer(abeta)
    bbeta <- as.integer(bbeta)
    ptop <- as.double(rep(0, NTOP))
    sigtop <- as.double(rep(0, NTOP))
    nftop <- as.integer(rep(0, NTOP))
    jtop <- matrix(0, nrow = NTOP, ncol = MXFAC)
    dimnames(jtop) <- list(paste("M", seq(NTOP), sep = ""), paste("x", 
        seq(MXFAC), sep = ""))
    storage.mode(jtop) <- "integer"
    prob <- as.double(rep(0, (COLS + 1)))
    names(prob) <- c("none", faclab)
    ind <- as.integer(-1)
    lst <- .Fortran("obm", X, Y, N, COLS, abeta, bbeta, BLKS, MXFAC, MXINT, 
        NTOP, mdcnt, ptop, nftop, jtop, prob, sigtop, ind , PACKAGE = "OBsMD")
    names(lst) <- c("X", "Y", "N", "COLS", "abeta", "bbeta", "BLKS", "MXFAC", "MXINT", 
        "NTOP", "mdcnt", "ptop",  "nftop", "jtop", "prob", 
        "sigtop","ind")
    invisible(structure(lst, class = c("OBsProb", class(lst))))
}
