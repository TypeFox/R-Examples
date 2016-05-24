MD <-
function (X, y, nFac, nBlk = 0, mInt = 3, g = 2, nMod, p, s2, 
    nf, facs, nFDes = 4, Xcand, mIter = 20, nStart = 5, startDes = NULL, 
    top = 20, eps = 1e-05) 
{
    if (nFac + nBlk != ncol(X)) 
        stop("nFac + nBlk != ncol(X)")
    if (nFac + nBlk != ncol(Xcand)) 
        stop("nFac + nBlk != ncol(Xcand)")
    if (ncol(Xcand) != ncol(X)) 
        stop("ncol(Xcand) != ncol(X)")
    ITMAX <- as.integer(mIter)
    N0 <- as.integer(nrow(X))
    NRUNS <- as.integer(nFDes)
    N <- as.integer(nrow(Xcand))
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    Y <- as.double(y)
    GAMMA <- as.double(g[1])
    GAM2 <- as.double(0)
    if (length(g) > 1) 
        GAM2 <- as.double(g[2])
    COLS <- as.integer(nFac)
    BL <- as.integer(nBlk)
    CUT <- as.integer(mInt)
    GAMMA <- as.double(g[1])
    if (length(g) == 1) {
        IND <- as.integer(0)
    }
    else {
        IND <- as.integer(1)
        GAM2 <- as.double(g[2])
    }
    Xcand <- as.matrix(Xcand)
    storage.mode(Xcand) <- "double"
    NM <- as.integer(nMod)
    P <- as.double(as.numeric(p))
    SIGMA2 <- as.double(as.numeric(s2))
    NF <- as.integer(as.numeric(nf))
    MNF <- as.integer(max(NF))
    JFAC <- as.matrix(facs)
    storage.mode(JFAC) <- "integer"
    if (is.null(startDes)) {
        if (is.null(nStart)) 
            stop("nStart needed when startDes is NULL")
        INITDES <- as.integer(1)
        NSTART <- as.integer(nStart)
        MBEST <- matrix(0, nrow = NSTART, ncol = NRUNS)
        storage.mode(MBEST) <- "integer"
    }
    else {
        INITDES <- as.integer(0)
        startDes <- as.matrix(startDes)
        NSTART <- as.integer(nrow(startDes))
        if (ncol(startDes) != NRUNS) 
            stop("ncol(startDes) should be nFDes")
        MBEST <- as.matrix(startDes)
        storage.mode(MBEST) <- "integer"
    }
    NTOP <- as.integer(top)
    TOPD <- as.double(rep(0, NTOP))
    TOPDES <- matrix(0, nrow = NTOP, ncol = NRUNS)
    dimnames(TOPDES) <- list(seq(top), paste("r", seq(NRUNS), 
        sep = ""))
    storage.mode(TOPDES) <- "integer"
    EPS <- as.double(eps)
    flag <- as.integer(-1)
    lst <- .Fortran("md", NSTART, NRUNS, ITMAX, INITDES, N0, 
        IND, X, Y, GAMMA, GAM2, BL, COLS, N, Xcand, NM, P, SIGMA2, 
        NF, MNF, JFAC, CUT, MBEST, NTOP, TOPD, TOPDES, EPS, flag, 
        PACKAGE = "BsMD")
    names(lst) <- c("NSTART", "NRUNS", "ITMAX", "INITDES", "N0", 
        "IND", "X", "Y", "GAMMA", "GAM2", "BL", "COLS", "N", 
        "Xcand", "NM", "P", "SIGMA2", "NF", "MNF", "JFAC", "CUT", 
        "MBEST", "NTOP", "TOPD", "TOPDES", "EPS", "flag")
    invisible(structure(lst, class = c("MD", class(lst))))
}
