GAopt <- function (OF, algo = list(), ...) {
    algoD <- list(nB = NA,     ### bits per solution
                  nP = 50L,    ### population size
                  nG = 300L,   ### number of generations
                  prob = 0.01, ### probability of mutation
                  pen = NULL,
                  repair = NULL,
                  loopOF = TRUE,
                  loopPen = TRUE,
                  loopRepair = TRUE,
                  methodOF = "loop",    ## 'vectorised','snow','multicore'
                  cl = NULL,            ##
                  mc.control = list(),  ##
                  printDetail = TRUE,
                  printBar = TRUE,
                  initP = NULL,
                  storeF = TRUE,
                  storeSolutions = FALSE,
                  crossover = c("onePoint", "uniform")
                  )

    checkList(algo, algoD)
    algoD[names(algo)] <- algo

    if (!exists(".Random.seed", envir = .GlobalEnv,
                inherits = FALSE))
        state <- NA else state <- .Random.seed

    methodOF <- tolower(algoD$methodOF)
    if (methodOF == "vectorised" && identical(algoD$loopOF, TRUE))
        algoD$loopOF <- FALSE

    dosnow <- FALSE
    domc <- FALSE
    cl <- algoD$cl

    if (methodOF == "snow") {
        if (is.null(cl)) {
            method <- "loop"
            warning("no cluster ", sQuote("cl"),
                    " passed for method ", sQuote("snow"),
                    ": will use method ", sQuote("loop"))
        } else {
            dosnow <- TRUE
            if (is.numeric(cl)) {
                cl <- makeCluster(c(rep("localhost", cl)), type = "SOCK")
                on.exit(stopCluster(cl))
            }
        }
    } else if (methodOF == "multicore") {
        domc <- TRUE
        mc.settings <- mcList(algoD$mc.control)
    }

    printDetail <- algoD$printDetail
    printBar <- algoD$printBar
    if (printBar && printDetail > 1)
        printBar <- FALSE

    if (!is.function(OF))
        stop("'OF' must be a function")
    OF1 <- function(x)
        OF(x, ...)
    Pe1 <- function(x)
        algoD$pen(x, ...)
    Re1 <- function(x)
        algoD$repair(x, ...)


    if (is.null(algoD$nB))
        stop("'algo$nB' must be specified")
    if (algoD$prob > 1 || algoD$prob < 0)
        stop("'algo$prob' must be between 0 and 1")

    nB <- makeInteger(algoD$nB, 'algo$nB')
    nP <- makeInteger(algoD$nP, 'algo$nP')
    nG <- makeInteger(algoD$nG, 'algo$nG')
    lP <- vector(mode = "list", length = nP)

    crossover <- tolower(algoD$crossover[1L])
    crossOver1 <- FALSE
    crossOver2 <- FALSE
    if (crossover == "onepoint") {
        crossOver <- function(x,y) {
            cutoff <- sample.int(nB-1L, 1L) + 1L
            ii <- cutoff:nB
            x[ii] <- y[ii]
            x
        }
        crossOver1 <- TRUE
    } else if (crossover == "uniform") {
        crossOver <- function(x,y) {
            ii <- runif(nB) > 0.5
            x[ii] <- y[ii]
            x
        }
        crossOver2 <- TRUE
    } else {
        stop("unknown crossover type")
    }

    snP <- seq_len(nP)
    snP1 <- c(nP, snP[-nP])

    vF <- numeric(nP)
    vF[] <- NA
    vFc <- vP <- vF

    if (algoD$storeF)
        Fmat <- array(NA, c(nG, nP)) else Fmat <- NA

    ## create population
    if (is.null(algoD$initP)) {
        mP <- array(sample.int(2L, nB * nP, replace = TRUE) - 1L,
                    dim = c(nB,nP))
        storage.mode(mP) <- "logical"
    } else {
        if (is.function(algoD$initP))
            mP <- algoD$initP()
        else
            mP <- algoD$initP
        if (mode(mP) != "logical") {
            storage.mode(mP) <- "logical"
            warning("'initP' is not of mode logical; 'storage.mode(initP)' will be tried")
        }
    }
    if (algoD$storeSolutions)
        xlist <- list(P = vector("list", length = nG), initP = mP) else xlist <- NA

    ## repair initial population
    if (!is.null(algoD$repair)) {
        if (algoD$loopRepair) {
            for (s in snP) mP[, s] <- Re1(mP[, s])
        } else {
            mP <- Re1(mP)
        }
    }

    ## evaluate initial population
    if (algoD$loopOF) {
        if (dosnow) {
            for (s in snP)
                lP[[s]] <- mP[ ,s]
            vF <- unlist(clusterApply(cl, lP, OF, ...))
        } else if (domc) {
            for (s in snP)
                lP[[s]] <- mP[ ,s]
            vF <- unlist(mclapply(lP, OF, ...,
                                  mc.preschedule = mc.settings$mc.preschedule,
                                  mc.set.seed = mc.settings$mc.set.seed,
                                  mc.silent = mc.settings$mc.silent,
                                  mc.cores = mc.settings$mc.cores,
                                  mc.cleanup = mc.settings$mc.cleanup))
        } else {
            for (s in snP)
                vF[s] <- OF1(mP[, s])
        }
    } else {
        vF <- OF1(mP)
    }
    if (!is.null(algoD$pen)) {
        if (algoD$loopPen) {
            for (s in snP) vP[s] <- Pe1(mP[, s])
        } else {
            vP <- Pe1(mP)
        }
        vF <- vF + vP
    }
    if (printBar)
        whatGen <- txtProgressBar(min = 1, max = nG, style = 3)
    if (printDetail)
        cat("\nGenetic Algorithm.\n")
    for (g in seq_len(nG)) {
        if (printBar)
            setTxtProgressBar(whatGen, value = g)

        ## create children ...
        mC <- array(NA, dim = dim(mP))

        ## ... through crossover
        o <- sample.int(nP)
        oo <- o[snP1]
        if (crossOver1) {
            for (s in snP)
                mC[ ,s] <- crossOver(mP[ ,o[s]],mP[ ,oo[s]])
        } else if (crossOver2) {
            mC <- crossOver(mP[ ,o],mP[ ,oo])
        }

        ## ... through mutation
        mutate <- runif(nB*nP) < algoD$prob
        mC[mutate] <- !mC[mutate]

        if (!is.null(algoD$repair)) {
            if (algoD$loopRepair) {
                for (s in snP) mC[ ,s] <- Re1(mC[ ,s])
            } else {
                mC <- Re1(mC)
            }
        }
        if (algoD$loopOF) {
            if (dosnow) {
                for (s in snP) lP[[s]] <- mC[ ,s]
                vFc <- unlist(clusterApply(cl, lP, OF, ...))
            } else if (domc) {
                for (s in snP) lP[[s]] <- mC[ ,s]
                vFc <- unlist(mclapply(lP, OF, ...,
                                    mc.preschedule = mc.settings$mc.preschedule,
                                    mc.set.seed = mc.settings$mc.set.seed,
                                    mc.silent = mc.settings$mc.silent,
                                    mc.cores = mc.settings$mc.cores,
                                    mc.cleanup = mc.settings$mc.cleanup
                                    )
                              )
            } else {
                for (s in snP) vFc[s] <- OF1(mC[, s])
            }
        } else {
            vFc <- OF1(mC)
        }
        if (!is.null(algoD$pen)) {
            if (algoD$loopPen) {
                for (s in snP) vP[s] <- Pe1(mC[ ,s])
            } else {
                vP <- Pe1(mC)
            }
            vFc <- vFc + vP
        }
        ## pairwise comparison
        logik <- vFc < vF
        mP[, logik] <- mC[, logik]
        vF[logik] <- vFc[logik]
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
    }
    if (printBar)
        close(whatGen)
    sGbest <- min(vF)
    sgbest <- which.min(vF)[1L]

    if (printDetail)
        cat("Best solution has objective function value ",
            prettyNum(sGbest), " ;",
            "\nstandard deviation of OF in final population is ",
            prettyNum(sd(vF)), " .\n\n", sep = "")

    list(xbest = mP[, sgbest], OFvalue = sGbest, popF = vF,
         Fmat = Fmat, xlist = xlist, initial.state = state)
}
