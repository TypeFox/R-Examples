test.TAopt <- function() {

    ## TA should come close to the minimum
    size <- 5L
    x0 <- runif(size)
    xTRUE <- runif(size)
    data <- list(xTRUE = xTRUE,
                 step = 0.02)
    OF <- function(x, data)
        max(abs(x - data$xTRUE))
    neighbour <- function(x, data)
        x + runif(length(data$xTRUE))*data$step - data$step/2

    algo <- list(q = 0.05, nS = 1000L, nT = 10L,
                 neighbour = neighbour, x0 = x0,
                 printBar = FALSE,
                 printDetail = FALSE,
                 storeSolutions = TRUE,
                 storeF = TRUE)
    res <- TAopt(OF, algo = algo, data = data)
    checkTrue(res$OFvalue < 0.005)

    ## check 'xlist'
    checkTrue(length(res$xlist[[1L]])==length(res$xlist[[2L]]))
    checkEquals(dim(do.call(rbind, res$xlist[[1L]])),
                dim(do.call(rbind, res$xlist[[2L]])))
    checkEquals(dim(res$Fmat), c(algo$nS * algo$nT, 2L))

    ## check 'Fmat': xn and xc must not all be the same
    checkTrue(!isTRUE(all.equal(res$Fmat[ ,1L],res$Fmat[ ,2L])))

    tmp <- res$Fmat[ ,1L]==res$Fmat[ ,2L]
    checkEquals(res$Fmat[tmp ,1L],
                apply(do.call(rbind, res$xlist[[1L]])[tmp, ], 1,OF, data))

    ## length(returned thresholds) == specified length(thresholds)
    checkTrue(length(res$vT) == algo$nT)

    ## specified thresholds are used
    algo$vT <- c(0.1,0.05,0)
    algo$nS <- 10L
    res <- TAopt(OF, algo = algo, data = data)
    checkEqualsNumeric(res$vT,algo$vT)

    ## stepUp is used
    algo$stepUp <- 2L
    res <- TAopt(OF, algo = algo, data = data)
    checkEqualsNumeric(res$vT, rep(algo$vT, 3L))

    ## scale is used
    algo$stepUp <- 0L
    algo$scale <- 1.5
    res <- TAopt(OF, algo = algo, data = data)
    checkEqualsNumeric(res$vT, algo$scale*c(0.1,0.05,0))

    ## q is zero
    algo <- list(q = 0, nS = 10L, nT = 15L,
                 neighbour = neighbour, x0 = x0,
                 printBar = FALSE,
                 printDetail = FALSE)
    res <- TAopt(OF, algo = algo, data = data)
    checkEqualsNumeric(res$vT, numeric(algo$nT))
    checkEquals(length(res$vT), algo$nT)

    ## check printDetail
    algo <- list(q = 0, nS = 10L, nT = 5L,
                 neighbour = neighbour, x0 = x0,
                 printBar = TRUE,
                 printDetail = 5)
    res <- capture.output(ignore <- TAopt(OF, algo = algo, data = data))
    checkEquals(sum(grepl("Best solution", res)), 11L)
}
