BsProb <-
function (X, y, blk = 0, mFac = 3, mInt = 2, p = 0.25, g = 2, 
    ng = 1, nMod = 10) 
{
    X <- as.matrix(X)
    y <- unlist(y)
    if (length(y) != nrow(X)) 
        stop("X and y should have same number of observations")
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
    PI <- as.double(p)
    if (length(g) == 1) {
        INDGAM <- as.integer(0)
        GAMMA <- as.double(g)
        NGAM <- as.integer(1)
        INDG2 <- as.integer(0)
        GAM2 <- as.double(0)
    }
    else {
        if (ng == 1) {
            INDGAM <- as.integer(0)
            GAMMA <- as.double(c(g[1], g[2]))
            NGAM <- as.integer(1)
            INDG2 <- as.integer(1)
            GAM2 <- as.double(g[2])
        }
        else {
            INDGAM <- as.integer(1)
            GAMMA <- as.double(seq(min(g), max(g), length = ng))
            NGAM <- as.integer(ng)
            INDG2 <- as.integer(0)
            GAM2 <- as.double(0)
        }
    }
    NTOP <- as.integer(nMod)
    mdcnt <- as.integer(0)
    ptop <- as.double(rep(0, NTOP))
    sigtop <- as.double(rep(0, NTOP))
    nftop <- as.integer(rep(0, NTOP))
    jtop <- matrix(0, nrow = NTOP, ncol = MXFAC)
    dimnames(jtop) <- list(paste("M", seq(NTOP), sep = ""), paste("x", 
        seq(MXFAC), sep = ""))
    storage.mode(jtop) <- "integer"
    del <- as.double(0)
    sprob <- as.double(rep(0, (COLS + 1)))
    names(sprob) <- c("none", faclab)
    pgam <- as.double(rep(0, NGAM))
    prob <- matrix(0, nrow = (1 + COLS), ncol = NGAM)
    dimnames(prob) <- list(c("none", paste("x", 1:COLS, sep = "")), 
        seq(NGAM))
    storage.mode(prob) <- "double"
    ind <- as.integer(-1)
    lst <- .Fortran("bm", X, Y, N, COLS, BLKS, MXFAC, MXINT, 
        PI, INDGAM, INDG2, GAM2, NGAM, GAMMA, NTOP, mdcnt, ptop, 
        sigtop, nftop, jtop, del, sprob, pgam, prob, ind, PACKAGE = "BsMD")
    names(lst) <- c("X", "Y", "N", "COLS", "BLKS", "MXFAC", "MXINT", 
        "PI", "INDGAM", "INDG2", "GAM2", "NGAM", "GAMMA", "NTOP", 
        "mdcnt", "ptop", "sigtop", "nftop", "jtop", "del", "sprob", 
        "pgam", "prob", "ind")
    invisible(structure(lst, class = c("BsProb", class(lst))))
}
