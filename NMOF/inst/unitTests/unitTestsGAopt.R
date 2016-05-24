## -*- truncate-lines: t; -*-
test.GAopt <- function() {

    size <- 20L; y <- runif(size) > 0.5; x <- runif(size) > 0.5
    data <- list(y = y); rm(y)
    algo <- list(nB = size, nP = 25L, nG = 150L, prob = 0.002,
                 printBar = FALSE, printDetail = FALSE,
                 crossover = "uniform",
                 storeSolutions = TRUE)
    OF <- function(x, data) sum(x != data$y)
    solGA <- GAopt(OF, algo = algo, data = data)
    ## should solve problem (solution is zero)
    checkEquals(solGA$OFvalue, 0)
    ## solution length == size
    checkEquals(length(solGA$xbest), size)
    ## population size
    checkEquals(length(solGA$popF), algo$nP)
    ## one matrix per generation
    checkEquals(length(solGA$xlist[[1L]]), algo$nG)
    ## each matrix in xlist should have dim c(size, algo$nP)
    checkEquals(sapply(solGA$xlist[[1L]], dim),
                array(1, dim = c(2L,algo$nG)) * c(size, algo$nP))
    ## there is no third element in 'xlist'
    checkException(length(solGA$xlist[[3L]]), silent = TRUE)

    ## error -- wrong size of initP
    algo$initP <- array(FALSE, dim = c(20,10))
    checkException(res <- GAopt(OF = OF, algo), silent = TRUE)
    algo$initP <- function() array(FALSE, dim = c(20L,20L))
    checkException(res <- GAopt(OF = OF, algo), silent = TRUE)
    algo$initP <- NULL

    ## error -- prob > 1, prob < 0
    algo$prob <- 2
    checkException(res <- GAopt(OF = OF, algo), silent = TRUE)
    algo$prob <- -0.1
    checkException(res <- GAopt(OF = OF, algo), silent = TRUE)

    ## error -- wrong crossover type
    algo <- list(nB = size, nP = 25L, nG = 150L, prob = 0.002,
                 printBar = FALSE, printDetail = FALSE,
                 crossover = "twoPoint")
    checkException(solGA <- GAopt(OF, algo = algo, data = data),
                   silent = TRUE)

    ## error -- OF not a function
    algo <- list(nB = size, nP = 25L, nG = 150L, prob = 0.002,
                 printBar = FALSE, printDetail = FALSE,
                 crossover = "uniform",
                 storeSolutions = TRUE)
    a <- numeric(10L)
    checkException(solGA <- GAopt(a, algo = algo, data = data),
                   silent = TRUE)

    ## check repair
    data$resample <- function(x, ...) x[sample.int(length(x), ...)]
    repairK <- function(x, data) {
        sx <- sum(x)
        if (sx > data$kmax) {
            i <- data$resample(which(x), sx - data$kmax)
            x[i] <- FALSE
        }
        x
    }
    repairK2 <- function(x, data) {
        sx <- colSums(x)
        whichCols <- which(sx > data$kmax)
        for (j in seq(along.with = whichCols)) {
            jj <- whichCols[j]
            i <- data$resample(which(x[ ,jj]), sx[jj] - data$kmax)
            x[i, jj] <- FALSE
        }
        x
    }
    data$kmax <- 5
    tempP <- array(TRUE, dim = c(20,10))
    checkTrue(all(colSums(repairK2(tempP,data))<=data$kmax))

    algo$repair <- repairK
    solGA <- GAopt(OF, algo = algo, data = data)
    checkTrue(sum(solGA$xbest)<=data$kmax)

    algo$repair <- repairK2; algo$loopRepair <- FALSE
    solGA <- GAopt(OF, algo = algo, data = data)
    checkTrue(sum(solGA$xbest)<=data$kmax)
    ## redo tests: repairing should not damage anything
    ## solution length == size
    checkEquals(length(solGA$xbest), size)
    ## population size
    checkEquals(length(solGA$popF), algo$nP)
    ## one matrix per generation
    checkEquals(length(solGA$xlist[[1L]]), algo$nG)
    ## each matrix in xlist should have dim c(size, algo$nP)
    checkEquals(sapply(solGA$xlist[[1L]], dim),
                array(1, dim = c(2L,algo$nG)) * c(size, algo$nP))
    ## there is no third element in 'xlist'
    checkException(length(solGA$xlist[[3L]]), silent = TRUE)

    ## check: do _not_ store solutions
    size <- 10L; OF <- function(x, y) sum(x != y)
    y <- runif(size) > 0.5 ## the true solution
    algo <- list(nB = size, nP = 20L, nG = 10L, prob = 0.002,
                 printBar = FALSE, printDetail = FALSE)
    xlist <- GAopt(OF, algo = algo, y = y)$xlist
    checkTrue(is.na(xlist))
    checkEquals(length(xlist), 1L)

    ## warning changed to error/reset on.exit
    ## checks 'snow': 'cl' not supplied
    ## op <- options("warn"); on.exit(options(op))
    ## options(warn = 2)
    ## size <- 10L
    ## OF <- function(x, y) sum(x != y)
    ## y <- runif(size) > 0.5 ## the true solution
    ## algo <- list(nB = size, nP = 20L, nG = 10L, prob = 0.002,
    ##              printBar = FALSE, methodOF = "snow")
    ## if (require("snow"))
    ##     checkException(sol <- GAopt(OF, algo = algo, y = y), silent = TRUE)
}
