
context("est_coi_um")

test_that("estimated intensity is correct in simple cases", {

    xoloc <- list(c(1, 2, 3), 2, c(1, 3))
    sclength <- c(4, 4, 4)
    centromere <- c(2, 2, 2)

    # estimated intensity with window = 0.05
    z05 <- est.coi.um(xoloc, sclength, centromere, intwindow=0.05, intloc=seq(0, 1, len=500))
    pos05 <- z05$intensity[,1]
    expected_intensity <- rep(0, length(pos05))
    expected_intensity[abs(pos05 - 0.25) <= 0.05/2 | abs(pos05 - 0.5) <= 0.05/2 | abs(pos05 - 0.75) <= 0.05/2] <- (2/3)/0.05
    expect_equal(z05$intensity[,2], expected_intensity, tolerance=1e-12)

    # estimated intensity with window = 0.10
    z10 <- est.coi.um(xoloc, sclength, centromere, intwindow=0.10, intloc=seq(0, 1, len=500))
    pos10 <- z10$intensity[,1]
    expected_intensity <- rep(0, length(pos10))
    expected_intensity[abs(pos10 - 0.25) <= 0.10/2 | abs(pos10 - 0.5) <= 0.10/2 | abs(pos10 - 0.75) <= 0.10/2] <- (2/3)/0.10
    expect_equal(z10$intensity[,2], expected_intensity, tolerance=1e-12)
})

test_that("same thing, with one cell having no crossovers", {

    xoloc <- list(c(1, 2, 3), 2, numeric(0), c(1, 3))
    sclength <- c(4, 4, 4, 4)
    centromere <- c(2, 2, 2, 2)

    # estimated intensity with window = 0.05
    z05 <- est.coi.um(xoloc, sclength, centromere, intwindow=0.05, intloc=seq(0, 1, len=500))
    pos05 <- z05$intensity[,1]
    expected_intensity <- rep(0, length(pos05))
    expected_intensity[abs(pos05 - 0.25) <= 0.05/2 | abs(pos05 - 0.5) <= 0.05/2 | abs(pos05 - 0.75) <= 0.05/2] <- 0.5/0.05
    expect_equal(z05$intensity[,2], expected_intensity, tolerance=1e-12)

    # estimated intensity with window = 0.10
    z10 <- est.coi.um(xoloc, sclength, centromere, intwindow=0.10, intloc=seq(0, 1, len=500))
    pos10 <- z10$intensity[,1]
    expected_intensity <- rep(0, length(pos10))
    expected_intensity[abs(pos10 - 0.25) <= 0.10/2 | abs(pos10 - 0.5) <= 0.10/2 | abs(pos10 - 0.75) <= 0.10/2] <- 0.5/0.10
    expect_equal(z10$intensity[,2], expected_intensity, tolerance=1e-12)
})

test_that("same thing, that numeric(0) being a NULL", {

    xoloc <- list(c(1, 2, 3), 2, NULL, c(1, 3))
    sclength <- c(4, 4, 4, 4)
    centromere <- c(2, 2, 2, 2)

    # estimated intensity with window = 0.05
    z05 <- est.coi.um(xoloc, sclength, centromere, intwindow=0.05, intloc=seq(0, 1, len=500))
    pos05 <- z05$intensity[,1]
    expected_intensity <- rep(0, length(pos05))
    expected_intensity[abs(pos05 - 0.25) <= 0.05/2 | abs(pos05 - 0.5) <= 0.05/2 | abs(pos05 - 0.75) <= 0.05/2] <- 0.5/0.05
    expect_equal(z05$intensity[,2], expected_intensity, tolerance=1e-12)

    # estimated intensity with window = 0.10
    z10 <- est.coi.um(xoloc, sclength, centromere, intwindow=0.10, intloc=seq(0, 1, len=500))
    pos10 <- z10$intensity[,1]
    expected_intensity <- rep(0, length(pos10))
    expected_intensity[abs(pos10 - 0.25) <= 0.10/2 | abs(pos10 - 0.5) <= 0.10/2 | abs(pos10 - 0.75) <= 0.10/2] <- 0.5/0.10
    expect_equal(z10$intensity[,2], expected_intensity, tolerance=1e-12)
})

test_that("estimated intensity with two groups and varying centromere pos and SC length", {

    xoloc <- list(c(1, 3, 4), numeric(0), 2, 4,
                  3, c(1,4), numeric(0), 2, numeric(0))
    sclength <- c(5, 1, 4, 5,
                  5, 5, 2, 4, 3)
    centromere <- c(2, 0.5, 2, 3,
                    2, 1.5, 1, 2, 1.5)
    group <- c("F", "F", "F", "F",
               "M", "M", "M", "M", "M")

    z <- est.coi.um(xoloc, sclength, centromere, group, intwindow=0.05, intloc=seq(0, 1, len=500))
    pos <- z$intensity[,1]

    # expected intensity for F group
    expintF <- rep(0, length(pos))
    loc <- c(0.25, 2/3, 5/6, 0.5, 0.75)
    for(x in loc)
        expintF[abs(pos - x) <= 0.05/2] <- (1/4)/0.05
    expect_equal(z$intensity[,2], expintF, tolerance=1e-12)

    # expected intensity for M group
    expintM <- rep(0, length(pos))
    loc <- c(2/3, 1/3, 6/7, 0.5)
    for(x in loc)
        expintM[abs(pos - x) <= 0.05/2] <- (1/5)/0.05
    expect_equal(z$intensity[,3], expintM, tolerance=1e-12)

})

