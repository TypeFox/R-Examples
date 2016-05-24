TAopt <- function(OF, algo = list(), ...) {
    algoD <- list(nD = 2000L, ## random steps for computing thresholds
                  nT = 10L,   ## number of thresholds
                  nS = 1000L, ## steps per threshold
                  q = 0.5,    ## starting quantile for thresholds
                  x0 = NULL,  ## initial solution
                  vT = NULL,  ## threshold sequence
                  neighbour = NULL,
                  printDetail = TRUE,
                  printBar = TRUE,
                  stepUp = 0L,
                  scale = 1,
                  storeF = TRUE,
                  storeSolutions = FALSE)

    checkList(algo, algoD)
    algoD[names(algo)] <- algo
    if (!exists(".Random.seed", envir = .GlobalEnv,
                inherits = FALSE))
        state <- NA else state <- .Random.seed


    ## user *must* specify the following
    if (is.null(algoD$neighbour))
        stop("specify a neighbourhood function 'algo$neighbour'")
    if (!is.function(algoD$neighbour))
        stop("'algo$neighbour' must be a function")
    if (is.null(algoD$x0))
        stop("specify start solution 'algo$x0'")
    if (is.function(algoD$x0)) ## evaluate x0 if function
        x0 <- algoD$x0() else x0 <- eval(algoD$x0)

    OF1 <- function(x)
        OF(x, ...)
    N1 <- function(x)
        algoD$neighbour(x, ...)

    printDetail <- algoD$printDetail
    printBar <- algoD$printBar
    if (printBar && printDetail > 1)
        printBar <- FALSE
    if (printDetail)
        cat("\nThreshold Accepting.\n")


    nT <- makeInteger(algoD$nT, "'algo$nT'")
    nS <- makeInteger(algoD$nS, "'algo$nS'")
    nD <- makeInteger(algoD$nD, "'algo$nD'")
    stepUp <- makeInteger(algoD$stepUp, "'algo$stepUp'", 0L)
    niter <- nS * nT * (stepUp+1L)

    ## compute thresholds
    if (is.null(algoD$vT)) {
        if (algoD$q < .Machine$double.eps^0.5) {
            vT <- numeric(nT)
        } else {
            if (printDetail) {
                cat("\nComputing thresholds ... ")
                flush.console()
                gc(FALSE)
                startTime <- proc.time()
            }
            if (printBar)
                whatGen <- txtProgressBar(min = 1, max = nD, style = 3,
                                          getOption("width")*0.9)
            xc  <- x0
            xcF <- OF1(xc)
            diffF <- numeric(nD)
            diffF[] <- NA
            for (i in seq_len(nD)){
                if (printBar)
                    setTxtProgressBar(whatGen, value = i)
                xn  <- N1(xc)
                xnF <- OF1(xn)
                diffF[i] <- abs(xcF - xnF)
                xc  <- xn
                xcF <- xnF
            }
            vT <- algoD$q * ((nT - 1L):0)/nT
            if (any(is.na(diffF)))
                stop("objective function evaluated to NA")
            vT <- quantile(diffF, vT, na.rm = FALSE)
            vT[nT] <- 0 ### set last threshold to zero
            if (printBar)
                close(whatGen)
            if (printDetail) {
                cat("OK.")
                endTime <- proc.time()
                cat("\nEstimated remaining running time:",
                    as.numeric(endTime[3L] - startTime[3L]) /
                    nD * niter * (stepUp + 1L),
                    "secs.\n\n")
                flush.console()
            }
        }
    } else {
        vT <- algoD$vT
    }
    if (stepUp > 0L)
        vT <- rep.int(vT, stepUp + 1L)
    if (algoD$scale < 0) {
        scale <- 0
        warning(sQuote("scale"), " set to 0 (",
                sQuote("scale"), " must be nonnegative)")
    } else {
        scale <- algoD$scale
    }
    vT <- vT * scale
    nT <- length(vT)
    niter <- nS * nT

    ## evaluate initial solution
    xc <- x0
    xcF <- OF1(xc)
    xbest <- xc
    xbestF <- xcF

    if (algoD$storeF) {
        Fmat <- array(NA, dim = c(niter, 2L))
        colnames(Fmat) <- c("xnF", "xcF")
    } else Fmat <- NA
    if (algoD$storeSolutions)
        xlist <- list(xn = vector("list", length = niter),
                      xc = vector("list", length = niter)) else xlist <- NA
    if (printDetail) {
        cat("\nRunning Threshold Accepting...\n")
        cat("Initial solution: ", prettyNum(xbestF),"\n")
        flush.console()
    }
    if (printBar)
        whatGen <- txtProgressBar(min = 1, max = niter, style = 3,
                                  getOption("width")*0.9)

    ## main algorithm
    counter <- 0L
    for (t in seq_len(nT)) {
        for (s in seq_len(nS)) {
            ## number of iterations
            counter <- counter + 1L

            xn <- N1(xc)
            xnF <- OF1(xn)
            if (xnF <= xcF + vT[t]) {
                xc <- xn
                xcF <- xnF
                if (xnF <= xbestF) {
                    xbest <- xn
                    xbestF <- xnF
                }
            }

            if (printBar)
                setTxtProgressBar(whatGen, value = counter)

            ## store OF values
            if (algoD$storeF) {
                Fmat[counter, 1L] <- xnF ## proposed sol.
                Fmat[counter, 2L] <- xcF ## accepted sol. (cummin=xbestF)
            }

            ## store solutions
            if (algoD$storeSolutions) {
                xlist[[c(1L, counter)]] <- xn
                xlist[[c(2L, counter)]] <- xc
            }

            ## print info
            if (printDetail > 1) {
                if (counter %% printDetail == 0L) {
                    cat("Best solution (iteration ", counter,
                        "/", niter, "): ",
                        prettyNum(xbestF),"\n", sep = "")
                    flush.console()
                }
            }

        }
    }
    if (printDetail)
        cat("Finished.\nBest solution overall: ",
            prettyNum(xbestF), "\n", sep = "")
    if (printBar)
        close(whatGen)

    ## return best solution
    list(xbest = xbest, OFvalue = xbestF,
         Fmat = Fmat, xlist = xlist, vT = vT,
         initial.state = state)
}

TA.info <- function(n = 0L) {
    e <- parent.frame(3L + n)
    step <- NA
    threshold <- NA
    iteration <- NA
    iteration.sampling <- NA
    if (exists("i", envir = e, inherits = FALSE))
        step <- get("i", envir = e, inherits = FALSE)
    if (exists("s", envir = e, inherits = FALSE))
        step <- get("s", envir = e, inherits = FALSE)
    if (exists("t", envir = e, inherits = FALSE))
        threshold <- get("t", envir = e, inherits = FALSE)
    if (exists("counter", envir = e, inherits = FALSE))
        iteration <- get("counter", envir = e, inherits = FALSE)
    list(iteration.sampling = iteration.sampling,
         iteration = iteration,
         step = step,
         threshold = threshold)
}


                                        # METHODS (not exported)

print.TAopt <- function(x, ...) {
    cat("Threshold Accepting.\n")
    cat(".. objective function value of solution: ", x$OFvalue, "\n")
}

plot.TAopt <- function(x, y, ...) {
    plot(x$vT, xlab = "Threshold", ylab = "Values")
    dev.new()
}


