OMD <-
function (X, y, nFac, nBlk = 0, mInt, nMod, optop, osigtop, onftop, ojtop, nFoll,
          Xcand, mIter, nStart, startDes, top = 20)
{
  #  if (nFac + nBlk != ncol(X)) 
  #      stop("nFac + nBlk != ncol(X)")
   # if (nFac + nBlk != ncol(Xcand)) 
    #    stop("nFac + nBlk != ncol(Xcand)")
    #if (ncol(Xcand) != ncol(X)) 
     #   stop("ncol(Xcand) != ncol(X)")
    ITMAX <- as.integer(mIter)
    N0 <- as.integer(nrow(X))
    NRUNS <- as.integer(nFoll)
    N <- as.integer(nrow(Xcand))
    X <- as.matrix(X)
    storage.mode(X) <- "double"
    Y <- as.double(y)
    COLS <- as.integer(nFac)
    BL <- as.integer(nBlk)
    CUT <- as.integer(mInt)
    Xcand <- as.matrix(Xcand)
    storage.mode(Xcand) <- "double"
    NM <- as.integer(nMod)
    P <- as.double(as.numeric(optop))
    SIGMA2 <- as.double(as.numeric(osigtop))
    NF <- as.integer(as.numeric(onftop))
    MNF <- as.integer(max(NF))
    JFAC <- as.matrix(ojtop)
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
    flag <- as.integer(-1)
    lst <- .Fortran("omd", NSTART, NRUNS, ITMAX, INITDES, N0, 
        X, Y, BL, COLS, N, Xcand, NM, P, SIGMA2, 
        NF, MNF, JFAC, CUT, MBEST, NTOP, TOPD, TOPDES, flag,   
        PACKAGE = "OBsMD")
    names(lst) <- c("NSTART", "NRUNS", "ITMAX", "INITDES", "N0", 
        "X", "Y", "BL", "COLS", "N", 
        "Xcand", "NM", "P", "SIGMA2", "NF", "MNF", "JFAC", "CUT", 
        "MBEST", "NTOP", "TOPD", "TOPDES", "flag")
    invisible(structure(lst, class = c("OMD", class(lst))))
}
