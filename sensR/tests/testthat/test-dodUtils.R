

## Load *ALL* objects in the sensR package to access non-exported
## functions etc.:
##   library(testthat)
##   library(devtools)
##   Path <- "/Users/rhbc/Documents/Rpackages/sensR/pkg/sensR"
##   load_all(Path)

context("Tests of dod utility functions")


##################################################################

test_that("dodSim returns data as expected", {
    ## Simulate data:
    set.seed(127)
    Data <- dodSim(d.prime=1, method="equi.prob")

    ## dput(Data)
    RES <- structure(c(25L, 18L, 22L, 22L, 33L, 30L, 20L, 30L),
                     .Dim = c(2L, 4L),
                     .Dimnames = list(c("same-pairs", "diff-pairs"),
                     c("1", "2", "3", "4")))
    expect_equal(Data, RES)

    ## Simulate data:
    set.seed(127)
    Data <- dodSim(d.prime=1, method="LR.max")
    expect_equivalent(Data, rbind(c(57, 30, 12, 1),
                                  c(51, 26, 17, 6)))
    ## Simulate data:
    set.seed(127)
    Data <- dodSim(d.prime=1, method="se.min")
    expect_equivalent(Data, rbind(c(51, 31, 15, 3),
                                  c(45, 26, 20, 9)))
    ## Simulate data:
    set.seed(127)
    Data <- dodSim(d.prime=1, method="user.defined", tau=1:3)
    expect_equivalent(Data, rbind(c(47, 38, 13, 2),
                                  c(42, 32, 19, 7)))
})

test_that("data can be numeric vectors or factors", {
    ## Simulate data:
    set.seed(127)
    data <- dodSim(d.prime=1, method="equi.prob")
    ## dput(readdoddata(data[1, ], data[2, ]))
    RES <- structure(list(same = c(25L, 22L, 33L, 20L),
                          diff = c(18L, 22L, 30L, 30L)),
                     .Names = c("same", "diff"))
    ## Data as vectors:
    expect_equal(readdoddata(data[1, ], data[2, ]), RES)
    ## Data as factor:
    same <- factor(rep.int(1:4, data[1, ]))
    diff <- factor(rep.int(1:4, data[2, ]))
    expect_equal(readdoddata(same, diff), RES)

    ## dod_fit:
    dp <- dod_fit(data[1, ], data[2, ])$d.prime
    dp2 <- dod_fit(same, diff)$d.prime
    expect_equal(dp, dp2)
    ## dod:
    dp <- dod(data[1, ], data[2, ])$d.prime
    dp2 <- dod(same, diff)$d.prime
    expect_equal(dp, dp2)
    ## dod_null_tau:
    tau <- dod_null_tau(data[1, ], data[2, ])
    tau2 <- dod_null_tau(same, diff)
    expect_equal(tau, tau2)
    ## dod_null:
    nll <- dod_null(data[1, ], data[2, ])
    nll2 <- dod_null(same, diff)
    expect_equal(nll, nll2)
    ## dod_nll:
    nll <- dod_nll(tau, dp, data[1, ], data[2, ])
    nll2 <- dod_nll(tau, dp, same, diff)
    expect_equal(nll, nll2)
})

test_that("nll functions return expected values", {
    ## Simulate data:
    set.seed(127)
    data <- dodSim(d.prime=1, method="equi.prob")
    same <- data[1, ]; diff <- data[2, ]
    ## dod_null_internal:
    ## dput(dod_null_internal(same, diff))
    RES <- 274.808886408307
    expect_equal(dod_null_internal(same, diff), RES)

    ## dput(dod_nll_internal(1:3, 1, same, diff))
    RES <- 346.038490942147
    expect_equal(dod_nll_internal(1:3, 1, same, diff), RES)
    expect_equal(dod_nll_all_internal(c(1:3, 1), same, diff), RES)
})

test_that("Integer warning can be turned off", {
    ## Simulate data:
    set.seed(127)
    data <- dodSim(d.prime=1, method="equi.prob")
    data <- data + runif(length(data))
    same <- data[1, ]; diff <- data[2, ]

    expect_warning(dod(same, diff), "non-integer counts")
    expect_true(givesWarnings(dod(same, diff)))

    ctrl <- dodControl(integer.tol=1)
    expect_false(givesWarnings(dod(same, diff, control=ctrl)))

    expect_warning(dod_fit(same, diff), "non-integer counts")
    expect_true(givesWarnings(dod_fit(same, diff)))

    ctrl <- dodControl(integer.tol=1)
    expect_false(givesWarnings(dod_fit(same, diff, control=ctrl)))
})

test_that("Warnings are counted correctly", {
    fun <- function() {
        warning("first warning")
        warning("second warning")
    }
    expect_false(givesWarnings(1:4))
    expect_true(givesWarnings(log(-1)))
    expect_true(givesWarnings(fun()))

    expect_true(countWarnings(1:4) == 0L)
    expect_true(countWarnings(log(-1)) == 1L)
    expect_true(countWarnings(fun()) == 2L)
})


##################################################################


