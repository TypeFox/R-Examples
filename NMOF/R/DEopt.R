DEopt <- function(OF, algo = list(), ...) {
    ## defaults for list 'algo'
    algoD <- list(nP = 50L,
                  nG = 300L,
                  F = 0.5, CR = 0.9,
                  min = NULL,
                  max = NULL,
                  pen = NULL,
                  repair = NULL,
                  loopOF = TRUE,
                  loopPen = TRUE,
                  loopRepair = TRUE,
                  methodOF = "loop",
                  printDetail = TRUE,
                  printBar = TRUE,
                  initP = NULL,
                  storeF = TRUE,
                  storeSolutions = FALSE,
                  minmaxConstr = FALSE)

    checkList(algo, algoD)
    algoD[names(algo)] <- algo
    if (!exists(".Random.seed", envir = .GlobalEnv,
                inherits = FALSE))
        state <- NA else state <- .Random.seed

    vmax <- as.vector(algoD$max)
    vmin <- as.vector(algoD$min)
    mm <- algoD$minmaxConstr
    if (is.null(vmin))
        stop("specify 'min' vector")
    if (is.null(vmax))
        stop("specify 'max' vector")
    if (length(vmax) != length(vmin))
        stop("'max' and 'min' have different lengths")
    if (!is.vector(vmax))
        stop("'max' must be a vector")
    if (!is.vector(vmin))
        stop("'min' must be a vector")
    if (any(vmin > vmax))
        stop("(at least) some 'max' < 'min'")

    d <- length(vmax)  ##  number of decision variables

    if ( algoD$CR > 1 || algoD$CR < 0 )
        stop("'CR' must be between 0 and 1")
    if ( any(algoD$F > 1) || any(algoD$F < 0) )
        warning("'F' is typically between 0 and 1")
    if ( length(algoD$F) > 1L && length(algoD$F) != d ) {
        warning("only first element of 'F' will be used")
        F <- algoD$F[1L]
    }
    F <- algoD$F[1L]

    printDetail <- algoD$printDetail
    printBar <- algoD$printBar

    if (printBar && printDetail > 1)
        printBar <- FALSE

    nG <- makeInteger(algoD$nG, "'algoD$nG'", 1L)
    nP <- makeInteger(algoD$nP, "'algoD$nP'", 1L)

    OF1 <- function(x)
        OF(x, ...)
    Pe1 <- function(x)
        algoD$pen(x, ...)
    Re1 <- function(x)
        algoD$repair(x, ...)

    snP <- seq_len(nP)
    snP1 <- c(nP, snP[-nP])

    ## MAIN ALGORITHM
    ## set up initial population
    vF <- numeric(nP); vF[] <- NA
    vPv <- vF; vFv <- vF
    if (algoD$storeF)
        Fmat <- array(NA, c(nG, nP)) else Fmat <- NA

    if (is.null(algoD$initP)) {
        mP <- vmin + mRU(d, nP) * (vmax - vmin)
    } else {
        if (is.function(algoD$initP))
            mP <- algoD$initP() else mP <- algoD$initP
        if (any(dim(mP) != c(d,nP)))
            stop("supplied population has wrong dimension")
    }
    if (algoD$storeSolutions)
        xlist <- list(P = vector("list", length = nG), initP = mP)
    else
        xlist <- NA

    ## evaluate initial population
    ## 1) repair (only required if initP not specified)
    if (mm && !is.null(algoD$initP))
        mP <- repair1c(mP, vmax, vmin)
    if (!is.null(algoD$repair)) {
        if (algoD$loopRepair){
            for (s in snP)
                mP[ ,s] <- Re1(mP[ ,s])
        } else {
            mP <- Re1(mP)
            if (!all(dim(mP) == c(d, nP)))
                stop("repair function returned population matrix of wrong size")
        }
    }
    ## 2) evaluate
    if (algoD$loopOF){
        for (s in snP)
            vF[s] <- OF1(mP[ ,s])
    } else {
        vF <- OF1(mP)
        if (length(vF) != nP)
            stop("objective function returned vector of wrong length")
    }
    ## 3) penalise
    if (!is.null(algoD$pen)) {
        if (algoD$loopPen) {
            for (s in snP) vPv[s] <- Pe1(mP[ ,s])
        } else {
            vPv <- Pe1(mP)
            if (length(vPv) != nP)
                stop("penalty function returned vector of wrong length")
        }
        vF <- vF + vPv
    }

    if (printBar)
        whatGen <- txtProgressBar (min = 1, max = nG, style = 3)
    if (printDetail)
        cat('\nDifferential Evolution.\n')

    for (g in seq_len(nG)) {
        if(printBar)
            setTxtProgressBar(whatGen, value = g)

        ## update population
        vI <- sample.int(nP)
        R1 <- vI[snP1]
        R2 <- R1[snP1]
        R3 <- R2[snP1]

        ## prelim. update
        mPv <- mP[ ,R1] + F * (mP[ ,R2] - mP[ ,R3])
        vI <- runif(d * nP) > algoD$CR
        mPv[vI] <- mP[vI]

        ## repair/evaluate/penalise updated population
        if (mm)
            mPv <- repair1c(mPv, vmax, vmin)
        if (!is.null(algoD$repair)) {
            if (algoD$loopRepair) {
                for (s in snP)
                    mPv[ ,s] <- Re1(mPv[ ,s])
            } else {
                mPv <- Re1(mPv)
            }
        }
        if (algoD$loopOF) {
            for (s in snP)
                vFv[s] <- OF1(mPv[ ,s])
        } else {
            vFv <- OF1(mPv)
        }
        if (!is.null(algoD$pen)) {
            if (algoD$loopPen){
                for (s in snP)
                    vPv[s] <- Pe1(mPv[ ,s])
            } else {
                vPv <- Pe1(mPv)
            }
            vFv <- vFv + vPv
        }
        ## find improvements
        logik <- vFv < vF
        mP[ ,logik] <- mPv[ ,logik]
        vF[logik] <- vFv[logik]
        if (algoD$storeF)
            Fmat[g, ] <- vF
        if (algoD$storeSolutions)
            xlist[[c(1L, g)]] <- mP

        ## print info
        if (printDetail > 1) {
            if (g %% printDetail == 0L) {
                cat("Best solution (iteration ", g, "/", nG, "): ",
                    prettyNum(min(vF)[1L]),"\n", sep = "")
                flush.console()
            }
        }

    } ## end of generations

    if (printBar)
        close(whatGen)

    ## return best solution
    sGbest <- min(vF); sgbest <- which.min(vF)[1L]

    if (printDetail)
        cat("Best solution has objective function value ",
            prettyNum(sGbest), " ;",
            "\nstandard deviation of OF in final population is ",
            prettyNum(sd(vF)), " .\n\n", sep = "")


    list(xbest = mP[ ,sgbest], OFvalue = sGbest,
         popF = vF, Fmat = Fmat, xlist = xlist,
         initial.state = state)
}

## if (0L) {
##     DEoptim <- function(OF, ..., lower = -Inf, upper = Inf,
##                         control = list()) {
##         stop("not operational yet")
##     }
##     PSoptim <- function(OF, ..., lower = -Inf, upper = Inf,
##                         control = list()) {
##         stop("not operational yet")
##     }
##     TAoptim <- function(OF, neighbour, ..., lower = -Inf, upper = Inf,
##                         control = list()) {
##         stop("not operational yet")
##     }
## }
