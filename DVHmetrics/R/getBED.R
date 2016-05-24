#####---------------------------------------------------------------------------
## implement recycling rule for function arguments
#####---------------------------------------------------------------------------

recycle <-
function(...) {
    dots <- list(...)
    maxL <- max(sapply(dots, length))
    lapply(dots, rep, length=maxL)
}

#####---------------------------------------------------------------------------
## BED linear quadratic model
#####---------------------------------------------------------------------------

getBED <-
function(D=NULL, fd=NULL, fn=NULL, ab=NULL) {
    UseMethod("getBED")
}

getBED.default <-
function(D=NULL, fd=NULL, fn=NULL, ab=NULL) {
    stopifnot(!is.null(fd), !is.null(ab))

    argL <- suppressWarnings(recycle(D, fd, fn, ab))
    D    <- argL[[1]]
    fd   <- argL[[2]]
    fn   <- argL[[3]]
    ab   <- argL[[4]]

    keepAB <- ab > 0
    keepFD <- fd > 0
    if(any(!keepAB)) { warning("'ab' must be > 0") }
    if(any(!keepFD)) { warning("'fd' must be > 0") }
    if(is.null(D)) {
        if(is.null(fn)) {
            stop("Either 'D' or 'fn' must be specified")
        } else {
            fn <- as.integer(fn)
            keepFN <- fn > 0
            if(any(!keepFN)) { warning("'fn' must be an integer > 0") }
            keep <- keepAB & keepFD & keepFN
            D    <- fn*fd
        }
    } else {
        if(!is.null(fn)) { warning("'fn' is ignored if 'D' is given") }
        keepD <- D > 0
        if(any(!keepD)) { warning("'D' must be > 0") }
        keep <- keepD & keepAB & keepFD
    }

    BED <- rep(NA_real_, times=length(D))
    BED[keep] <- D[keep] * (1 + (fd[keep]/ab[keep]))
    data.frame(BED=BED, fractDose=fd, ab=ab)
}

getBED.DVHs <-
function(D=NULL, fd=NULL, fn=NULL, ab=NULL) {
    stopifnot(!is.null(D), !is.null(fd))

    if(!is.null(fn)) { warning("'fn' is ignored for the 'DVHs' method") }
    if(length(fd) > 1L) {
        warning("Only first element of 'fd' will be used")
        fd <- fd[1]
    }

    if(length(ab) > 1L) {
        warning("Only first element of 'ab' will be used")
        ab <- ab[1]
    }

    D$dvh[ , "dose"] <- getBED(D=D$dvh[ , "dose"], fd=fd, ab=ab)$BED
    D
}

getBED.DVHLst <-
function(D=NULL, fd=NULL, fn=NULL, ab=NULL) {
    stopifnot(!is.null(D))

    BEDl <- Map(getBED, D, fd=list(fd), fn=list(fn), ab=list(ab))
    class(BEDl) <- "DVHLst"
    attr(BEDl, which="byPat") <- attributes(D)$byPat
    BEDl
}

getBED.DVHLstLst <-
function(D=NULL, fd=NULL, fn=NULL, ab=NULL) {
    stopifnot(!is.null(D))

    BEDll <- Map(getBED, D, fd=list(fd), fn=list(fn), ab=list(ab))
    class(BEDll) <- "DVHLstLst"
    attr(BEDll, which="byPat") <- attributes(D)$byPat
    BEDll
}
