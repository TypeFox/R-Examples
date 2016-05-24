#####---------------------------------------------------------------------------
## isoeffective dose linear quadratic model
#####---------------------------------------------------------------------------

getIsoEffD <-
function(D1=NULL, D2=NULL, fd1=NULL, fd2=NULL, ab=NULL) {
    UseMethod("getIsoEffD")
}

getIsoEffD.default <-
function(D1=NULL, D2=NULL, fd1=NULL, fd2=NULL, ab=NULL) {
    stopifnot(!is.null(D1), !is.null(fd1), !is.null(ab))

    argL <- suppressWarnings(recycle(D1, D2, fd1, fd2, ab))
    maxL <- max(sapply(argL, length))
    D1   <- argL[[1]]
    D2   <- argL[[2]]
    fd1  <- argL[[3]]
    fd2  <- argL[[4]]
    ab   <- argL[[5]]

    keepD1  <- D1  > 0
    keepFD1 <- fd1 > 0
    keepAB  <- ab  > 0

    if(any(!keepD1))  { warning("'D1' must be > 0")  }
    if(any(!keepFD1)) { warning("'fd1' must be > 0") }
    if(any(!keepAB))  { warning("'ab' must be > 0")  }

    if(is.null(D2)) {
        ## calculate D2 -> also need fd2 (special case: EQD2 with fd2=2)
        if(is.null(fd2)) { stop("'fd2' must be given to calculate 'D2'") }

        keepFD2 <- fd2 > 0
        if(any(!keepFD2)) { warning("'fd2' must be > 0") }
        
        keep <- keepD1 & keepFD1 & keepFD2 & keepAB
        D2   <- rep(NA_real_, times=maxL)
        D2[keep] <- D1[keep] * (fd1[keep] + ab[keep]) / (fd2[keep] + ab[keep])
        D2
    } else if(is.null(fd2)) {
        ## calculate fd2 -> also need D2
        if(is.null(D2)) { stop("'D2' must be given to calculate 'fd2'") }
    
        keepD2 <- D2 > 0
        if(any(!keepD2)) { warning("'D2' must be > 0") }

        keep <- keepD1 & keepD2 & keepFD1 & keepAB
        fd2  <- rep(NA_real_, times=maxL)
        fd2[keep] <- (D1[keep] / D2[keep]) * (fd1[keep] + ab[keep]) - ab[keep]
        fd2
    } else {
        warning("either 'D2' or 'fd2' must be NULL")
        rep(NA_real_, times=maxL)
    }
}

getIsoEffD.DVHs <-
function(D1=NULL, D2=NULL, fd1=NULL, fd2=NULL, ab=NULL) {
    stopifnot(!is.null(D1), !is.null(fd1), !is.null(fd2))
    if(!is.null(D2)) { warning("'D2' is ignored for the 'DVHs' method") }

    if(length(fd1) > 1L) {
        warning("Only first element of 'fd1' will be used")
        fd1 <- fd1[1]
    }

    if(length(fd2) > 1L) {
        warning("Only first element of 'fd2' will be used")
        fd2 <- fd2[1]
    }

    if(length(ab) > 1L) {
        warning("Only first element of 'ab' will be used")
        ab <- ab[1]
    }

    D1$dvh[ , "dose"] <- getIsoEffD(D1=D1$dvh[ , "dose"], fd1=fd1, fd2=fd2, ab=ab)
    D1
}

getIsoEffD.DVHLst <-
function(D1=NULL, D2=NULL, fd1=NULL, fd2=NULL, ab=NULL) {
    stopifnot(!is.null(D1))

    IsoEDl <- Map(getIsoEffD, D1, fd1=list(fd1), fd2=list(fd2), ab=list(ab))
    class(IsoEDl) <- "DVHLst"
    attr(IsoEDl, which="byPat") <- attributes(D1)$byPat
    IsoEDl
}

getIsoEffD.DVHLstLst <-
function(D1=NULL, D2=NULL, fd1=NULL, fd2=NULL, ab=NULL) {
    stopifnot(!is.null(D1))

    IsoEDll <- Map(getIsoEffD, D1, fd1=list(fd1), fd2=list(fd2), ab=list(ab))
    class(IsoEDll) <- "DVHLstLst"
    attr(IsoEDll, which="byPat") <- attributes(D1)$byPat
    IsoEDll
}
