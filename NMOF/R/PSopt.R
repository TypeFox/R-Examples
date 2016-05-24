PSopt <- function(OF, algo = list(), ...) {

    ## defaults
    algoD <- list(nP = 100L, nG = 500L,
                  c1 = 1, c2 = 1,
                  iner = 0.9, initV = 1, maxV = 1,
                  min = NULL, max = NULL,
                  pen = NULL, repair = NULL, changeV = NULL,
                  loopOF = TRUE, loopPen = TRUE,
                  loopRepair  = TRUE, loopChangeV = TRUE,
                  printDetail = TRUE, printBar = TRUE,
                  initP = NULL,
                  storeF = TRUE,
                  storeSolutions = FALSE,
                  minmaxConstr = FALSE
                  )

    checkList(algo, algoD)
    algoD[names(algo)] <- algo
    if (!exists(".Random.seed", envir = .GlobalEnv,
                inherits = FALSE))
        state <- NA else state <- .Random.seed


    vmax <- as.vector(algoD$max)
    vmin <- as.vector(algoD$min)
    mm <- algoD$minmaxConstr
    nP <- as.integer(algoD$nP)
    nG <- as.integer(algoD$nG)
    if(is.null(vmin))
        stop("specify 'min' vector")
    if(is.null(vmax))
        stop("specify 'max' vector")
    if(length(vmax) != length(vmin))
        stop("max/min have different lengths")
    if (!is.vector(vmax))
        stop("'max' must be a vector")
    if (!is.vector(vmin))
        stop("'min' must be a vector")
    if (any(vmin > vmax))
        stop("(at least) some max < min")
    printDetail <- algoD$printDetail
    printBar <- algoD$printBar
    if (printBar && printDetail > 1)
        printBar <- FALSE

    OF1 <- function(x)
        OF(x, ...)
    Pe1 <- function(x)
        algoD$pen(x, ...)
    Re1 <- function(x)
        algoD$repair(x, ...)
    cV1 <- function(x)
        algoD$changeV(x, ...)

    ## auxiliary functions
    pmax2 <- function(x1,x2)
        ( (x1 + x2) + abs(x1 - x2) ) / 2
    pmin2 <- function(x1,x2)
        ( (x1 + x2) - abs(x1 - x2) ) / 2

    if (algoD$storeF)
        Fmat <- array(NA, c(nG, nP)) else Fmat <- NA


    ## set up initial population and velocity
    d <- length(vmax)
    vF <- numeric(nP); vF[] <- NA; vPv <- vF
    if (is.null(algoD$initP)) {
        mP <- vmin + mRU(d, nP) * (vmax - vmin)
    } else {
        if (is.function(algoD$initP))
            mP <- algoD$initP() else mP <- algoD$initP
        if (any(dim(mP) != c(d, nP)))
            stop("supplied population has wrong dimension")

    }
    if (algoD$storeSolutions)
        xlist <- list(P     = vector("list", length = nG),
                      Pbest = vector("list", length = nG), initP = mP) else xlist <- NA
    mV <- algoD$initV * mRN(d,nP)

    ## evaluate initial population
    if (mm)
        mP <- repair1c(mP, vmax, vmin)
    if (!is.null(algoD$repair)) {
        if (algoD$loopRepair) {
            for(s in seq_len(nP))
                mP[ ,s] <- Re1(mP[ ,s])
        } else {
            mP <- Re1(mP)
            if (!all(dim(mP) == c(d, nP)))
                stop("repair function returned population matrix of wrong size")
        }
    }
    if (algoD$loopOF) {
        for(s in seq_len(nP)) vF[s] <- OF1(mP[,s])
    } else {
        vF <- OF1(mP)
        if (length(vF) != nP)
            stop("objective function returned vector of wrong length")

    }
    if (!is.null(algoD$pen)) {
        if(algoD$loopPen){
            for(s in seq_len(nP)) vPv[s] <- Pe1(mP[,s])
        } else {
            vPv <- Pe1(mP)
            if (length(vPv) != nP)
                stop("penalty function returned vector of wrong length")
        }
        vF <- vF + vPv
    }

    ## set up best solutions
    mPbest <- mP          ## matrix of 'personally best' solutions
    vFbest <- vF          ## vector of OF of best solutions
    sGbest <- min(vFbest) ## scalar: best OF-value
    sgbest <- which.min(vFbest)[1L] ## scalar: best solution (counter)

    ## start iterations
    if (printDetail)
        cat('\nParticle Swarm Optimisation.\n')
    if (printBar)
        whatGen <- txtProgressBar (min = 1, max = nG, style = 3)
    for (g in seq_len(nG)) {
        if(printBar)
            setTxtProgressBar(whatGen, value = g)
        ## update population
        mDV <- algoD$c1 * runif(d*nP) * (mPbest - mP) +
               algoD$c2 * runif(d*nP) * (mPbest[ ,sgbest] - mP)
        mV <- algoD$iner * mV + mDV
        mV <- pmin2(mV,  algoD$maxV)
        mV <- pmax2(mV, -algoD$maxV)
        if (!is.null(algoD$changeV)) {
            if (algoD$loopChangeV){
                for (s in seq_len(nP))
                    mV[ ,s] <- cV1(mV[ ,s])
            } else
            mV <- cV1(mV)
        }
        mP <- mP + mV

        ## evaluate updated population
        if (mm)
            mP <- repair1c(mP, vmax, vmin)
        if (!is.null(algoD$repair)) {
            if (algoD$loopRepair){
                for (s in seq_len(nP))
                    mP[ ,s] <- Re1(mP[ ,s])
            } else
            mP <- Re1(mP)
        }

        if (algoD$loopOF) {
            for (s in seq_len(nP)) vF[s] <- OF1(mP[ ,s])
        } else
        vF <- OF1(mP)

        if(!is.null(algoD$pen)) {
            if (algoD$loopPen){
                for (s in seq_len(nP)) vPv[s] <- Pe1(mP[ ,s])
            } else
            vPv <- Pe1(mP)
            vF <- vF + vPv
        }
        ## find improvements
        logik <- vF < vFbest
        mPbest[ ,logik] <- mP[ ,logik]
        vFbest[logik] <- vF[logik]
        ## find best solution
        if (min(vF) < sGbest) {
            sGbest <- min(vF)
            sgbest <- which.min(vF)[1L]
        }

        if (algoD$storeF)
            Fmat[g, ] <- vFbest
        if (algoD$storeSolutions) {
            xlist[[c(1L, g)]] <- mP
            xlist[[c(2L, g)]] <- mPbest
        }

        ## print info
        if (printDetail > 1) {
            if (g %% printDetail == 0L) {
                cat("Best solution (iteration ", g, "/", nG, "): ",
                    prettyNum(sGbest[1L]), "\n", sep = "")
                flush.console()
            }
        }

    } ## end generations

    if (printBar)
        close(whatGen)

    if (printDetail)
        cat("Best solution has objective function value ",
            prettyNum(sGbest), " ;",
            "\nstandard deviation of OF in final population is ",
            prettyNum(sd(vF)), " .\n\n", sep = "")

    ## return best solution
    list(xbest = mPbest[,sgbest], OFvalue = sGbest,
         popF = vFbest, Fmat = Fmat, xlist = xlist,
         initial.state = state)
}
