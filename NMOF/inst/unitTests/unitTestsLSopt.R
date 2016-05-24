test.LSopt <- function() {
    xTRUE <- runif(5)
    data <- list(xTRUE = xTRUE, step = 0.02)
    OF <- function(x, data) max(abs(x - data$xTRUE))
    neighbour <- function(x, data)
        x + runif(length(data$xTRUE))*data$step - data$step/2

    x0 <- runif(5)
    algo <- list(nS = 10000L,
                 neighbour = neighbour, x0 = x0,
                 printBar = FALSE, printDetail = FALSE,
                 storeF = TRUE, storeSolutions = TRUE)
    res <- LSopt(OF, algo = algo, data = data)
    checkTrue(res$OFvalue < 0.005)

    ## check 'xlist'
    checkTrue(length(res$xlist[[1L]])==length(res$xlist[[2L]]))
    checkEquals(dim(do.call(rbind, res$xlist[[1L]])),
                dim(do.call(rbind, res$xlist[[2L]])))
    checkEquals(dim(res$Fmat), c(algo$nS, 2L))

    ## check 'Fmat': xn and xc must not all be the same
    checkTrue(!isTRUE(all.equal(res$Fmat[ ,1L],res$Fmat[ ,2L])))

    tmp <- res$Fmat[ ,1L]==res$Fmat[ ,2L]
    checkEquals(res$Fmat[tmp ,1L],
                apply(do.call(rbind, res$xlist[[1L]])[tmp, ], 1,OF, data))

    ## check printDetail
    algo <- list(nS = 100L,
                 neighbour = neighbour, x0 = x0,
                 printBar = FALSE, printDetail = 10)
    res2 <- capture.output(res <- LSopt(OF, algo = algo, data = data))
    checkEquals(sum(grepl("Best solution", res2)), 11L)
}
