## -*- truncate-lines: t; -*-
## gridSearch
## bracketing
## xwGauss

test.gridSearch <- function() {
    testFun  <- function(x) x[1L] + x[2L]^2
    testFun1 <- function(x,k) x[1L] + x[2L]^2 + k
    testFun2 <- function(x) x[[1L]] + x[[2L]]^2
    testFun3 <- function(x) sum(x^2)
    testFun4 <- function(x) sum(sapply(x, `^`, 2L))

    ## use 1 -- specify all levels
    levels <- list(a = 1:2, b = 1:3)
    res <- gridSearch(testFun, levels, printDetail = FALSE)
    checkEquals(res$minfun,2)
    checkEquals(res$values,
                apply(as.matrix(expand.grid(levels)), 1L, testFun))

    ## use 2 -- specify lower, upper and npar
    lower <- 1; upper <- 3; npar <- 2
    res <- gridSearch(testFun, lower = lower, upper = upper,
                      npar = npar, printDetail = FALSE)
    checkEquals(res$minfun,2)

    ## use 3 -- specify lower, upper, npar and n
    lower <- 1; upper <- 3; npar <- 2; n <- 4
    res <- gridSearch(testFun, lower = lower, upper = upper,
                      npar = npar, n = n, printDetail = FALSE)
    checkEquals(res$minfun,2)

    ## use 4 -- specify lower, upper and n (same result as 'use 3')
    lower <- c(1,1); upper <- c(3,3); n <- 4
    res1 <- gridSearch(testFun, lower = lower, upper = upper, n = n)
    checkEquals(res, res1)

    ## use 5 -- specify lower, upper (auto-expanded) and n
    ## (same result as 'use 3')
    lower <- c(1,1); upper <- 3; n <- 4
    res2 <- gridSearch(testFun, lower = lower, upper = upper,
                       n = n, printDetail = FALSE)
    checkEquals(res1, res2)

    ## use 6 -- specify lower (auto-expanded), upper and n
    ## (same result as 'use 3')
    lower <- 1; upper <- c(3,3); n <- 4
    res3 <- gridSearch(testFun, lower = lower, upper = upper,
                       n = n, printDetail = FALSE)
    checkEquals(res2, res3)

    ## additional argument passed with '...'
    levels <- list(a = 1:2, b = 1:3); k <- 5
    res <- gridSearch(testFun1, levels, k = k, printDetail = FALSE)
    checkEquals(res$minfun, 7)

    ## 'keepNames'
    testFunNames  <- function(x) x[["a"]] + x[["b"]]^2
    levels <- list(a = 1:2, b = 1:3)
    res <- gridSearch(testFunNames, levels, printDetail = FALSE,
                      keepNames = TRUE)
    checkEquals(res$minfun, 2)
    res <- gridSearch(testFunNames, levels, printDetail = FALSE,
                      keepNames = TRUE, asList = TRUE)
    checkEquals(res$minfun, 2)
    testFunNames  <- function(x) x[["a"]] + x[["C"]]^2
    checkException(gridSearch(testFunNames, keepNames = TRUE,
                              levels, printDetail = FALSE),
                   silent = TRUE)

    ## ERROR -- length(upper) != length(lower)
    lower <- 1:3; upper <- 5:4; n <- 8
    checkException(gridSearch(fun = testFun,
                              lower = lower, upper = upper,
                              n = n, printDetail = FALSE),
                   silent = TRUE)

    ## ERROR -- upper < lower
    lower <- 1:3; upper <- 2; n <- 8
    checkException(gridSearch(fun = testFun,
                              lower = lower, upper = upper, n = n,
                              printDetail = FALSE),
                   silent = TRUE)

    ##
    lower <- 1:3; upper <- 5; n <- 8
    sol <- gridSearch(fun = testFun,
                      lower = lower, upper = upper,
                      n = n, printDetail = FALSE)
    checkEquals(sol$minlevels,1:3)

    ## compare levels with 'expand.grid' results
    levels <- list(a = 1:2, b = 1:3, c = 4:6)
    sol <- gridSearch(fun = testFun3, levels = levels, printDetail = FALSE)
    checkEquals(sol$minlevels,c(1,1,4))
    ## check levels
    l1 <- do.call("rbind", sol$levels)
    l2 <- as.matrix(expand.grid(levels))
    dimnames(l2) <- NULL
    checkEquals(l1,l2)

    ## error: must be greater than 1
    lower <- 1; upper <- 5; n <- 1
    sol <- checkException(gridSearch(fun = testFun3,
                                     lower = lower, upper = upper,
                                     n = n, printDetail = FALSE),
                          silent = TRUE)
    ##
    ## NSS fit
    tm <- seq_len(10); trueP <- c(2,1,1,5,1,3); y <- NSS(trueP, tm)
    qfit <- function(x, tm, y) {
        X <- NSSf(x[1], x[2], tm)
        mm <- lm(y ~ -1 + X)
        sum(abs(y - NSS(c(coef(mm), x[1], x[2]), tm)))
    }
    res <- gridSearch(qfit, y = y, tm = tm,
                      lower = 0.0, upper = 5, npar = 2L, n = 11L)
    ll <-  res$minlevels
    X <- NSSf(ll[1L], ll[2L], tm)
    model <- lm(y ~ -1 + X)
    checkEquals(sum(abs(y - NSS(c(coef(model), ll[1L], ll[2L]), tm))), 0)
}


## bracketing
test.bracketing <- function() {
    ## example ch. 11/p. 290
    res0 <- structure(c(0.3, 0.348, 0.444, 0.78,
                        0.324, 0.372, 0.468, 0.804), .Dim = c(4L, 2L))
    testFun <- function(x) cos(1/x^2)
    checkEquals(res0, bracketing(testFun, interval = c(0.3, 0.9), n = 26L))
    checkEquals(res0, bracketing(testFun, interval = c(0.3, 0.9), n = 26L),
                method = "vectorise")
    checkEquals(res0,bracketing(testFun, lower = 0.3, upper = 0.9, n = 26L))

    ## 'lower'/'upper' ignored if 'interval' is specified
    checkEquals(res0, bracketing(testFun, interval = c(0.3, 0.9),
                                lower = 0.1, upper = 0.99, n = 26L))
    checkEquals(bracketing(testFun, interval = c(0.1, 0.99)),
                bracketing(testFun, lower = 0.1, upper = 0.99))

    ## ERROR: lower < upper
    checkException(bracketing(testFun, lower = 0.3, upper = 0.3, n = 26L),
                   silent = TRUE)
    ## ERROR: no interval
    checkException(bracketing(testFun, n = 26L), silent = TRUE)
    checkException(bracketing(testFun, lower = 0.1, n = 26L), silent = TRUE)
    checkException(bracketing(testFun, upper = 0.1, n = 26L), silent = TRUE)

    ## no zero
    testFun <- function(x) 1
    res <- bracketing(testFun, interval = c(1,10), n = 10L)
    checkTrue(all(dim(res) == c(0L, 2L)))

    ## no zero either
    testFun <- function(x) 0
    res <- bracketing(testFun, interval = c(1,10), n = 10L)
    checkTrue(all(dim(res) == c(0L, 2L)))

    ## additional parameter
    testFun <- function(x,k) cos(1/x^2) + k
    res <- bracketing(testFun, k = 0, interval = c(0.3,0.9), n = 26L,
                      method = "vectorise")
    checkEquals(res0, res)
    res <- bracketing(testFun, k = 1, interval = c(0.3,0.9), n = 26L,
                      method = "vectorise")
    checkTrue(all(dim(res) == c(0L, 2L)))
}


## integration
test.xwGauss <- function() {
    ## http://dlmf.nist.gov/3.5

    tab <- list(nodes = 0, weights = 2)
    res <- xwGauss(n =  1L, "legendre")
    checkEquals(tab, res)

    tab <- structure(list(nodes =
                          c(-0.97390652851717, -0.865063366688983,
                            -0.679409568299024, -0.433395394129246,
                            -0.148874338981631, 0.148874338981632,
                            0.433395394129248, 0.679409568299025,
                            0.865063366688984, 0.973906528517172),
                          weights =
                          c(0.0666713443086876, 0.149451349150583,
                            0.219086362515984, 0.269266719309995,
                            0.29552422471475, 0.295524224714753,
                            0.269266719309995, 0.219086362515983,
                            0.149451349150581, 0.066671344308688)),
                     .Names = c("nodes", "weights"))
    res <- xwGauss(n =  10L, "legendre")
    checkEquals(tab, res)

    ##xwGauss(n = 200, "legendre")

    ##xwGauss(n =   1, "laguerre")
    ##xwGauss(n =  10, "laguerre")
    ##xwGauss(n = 200, "laguerre")

    ##xwGauss(n =   1, "hermite")
    ##xwGauss(n =  10, "hermite")
    ##xwGauss(n = 200, "hermite")

    checkException(xwGauss(n = 1, "yo"), silent = TRUE)
    checkException(xwGauss(n = 0, "legendre"), silent = TRUE)
}
