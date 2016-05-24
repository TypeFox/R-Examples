

## Load *ALL* objects in the sensR package to access non-exported
## functions etc.:
##   library(testthat)
##   library(devtools)
##   Path <- "/Users/rhbc/Documents/Rpackages/sensR/pkg/sensR"
##   load_all(Path)

context("Tests of the main dod function(s)")

test_that("We get an error when we expect one", {
    ## expect_true(FALSE) ## Should tricker and R CMD check error
})

test_that("dod estimates are stable", {
    ## Simulate some data and estimate:
    dprime <- 1
    set.seed(127)
    data <- dodSim(d.prime=dprime, sample.size=100, method="equi.prob")
    fm <- dod(same=data[1, ], diff=data[2, ])

    res <- unlist(fm[c("d.prime", "tau", "logLik", "p.value", "stat.value")])
    ## dput(res)
    RES <-
        structure(c(0.951548109556325, 0.429492707467134, 0.905080420721259,
                    1.805555643669, -273.418937317776, 0.0477271183728321, 1.6673026663032
                    ), .Names = c("d.prime", "tau1", "tau2", "tau3", "logLik", "p.value",
                       "stat.value"))

    expect_equal(res, RES)

    ## dput(coef(fm))
    coefExpect <-
        structure(c(0.951548109556325, 0.308805702582903, 0, 1.47386430866794),
                  .Dim = c(1L, 4L), .Dimnames = list("d.prime", c("Estimates",
                                    "Std. Error", "Lower", "Upper")))

    expect_equal(coef(fm), coefExpect)
})

test_that("dod works for extreme data settings", {
    ## Large d-prime:
    d <- 10
    set.seed(127)
    (data <- dodSim(d.prime=d, sample.size=100, method="equi"))
    expect_warning(dod(data[1, ], data[2, ]),
                   "Cannot assess convergence: non-finite gradient")
    expect_warning(dod_fit(data[1, ], data[2, ]),
                   "Cannot assess convergence: non-finite gradient")

    ## Almost complete separation:
    data <- rbind(c(55, 45, 2, 0),
                  c(0, 0, 48, 52))
    expect_warning(dod(data[1, ], data[2, ]),
                   "Cannot assess convergence: non-finite gradient")
    expect_warning(dod_fit(data[1, ], data[2, ]),
                   "Cannot assess convergence: non-finite gradient")
    data <- rbind(c(55, 45, 0, 1),
                  c(0, 0, 48, 52))
    expect_warning(dod(data[1, ], data[2, ]),
                   "Estimation failed with max")

    data <- rbind(c(55, 45, 0, 0),
                  c(0, 0, 48, 52))
    data[, 2:3] <- 0
    ##  ## This does not work, but probably it should:
    ##  expect_error(dod(data[1, ], data[2, ]),
    ##               "Complete data separation: DOD model not defined")
    ##  ## This does not work, but probably it should:
    ##  expect_error(dod_fit(data[1, ], data[2, ]),
    ##               "Complete data separation: DOD model not defined")

    expect_warning(dod(data[2, ], data[1, ]),
                   "d.prime < 0.01: standard errors are unavailable")
    ## This does not work, but probably it should:
    ##  expect_warning(dod_fit(data[2, ], data[1, ]),
    ##                 "d.prime < 0.01: standard errors are unavailable")

    ## Another example:
    d <- 5
    tau <- optimal_tau(d.prime=d, ncat=4, method="equi.prob")$tau
    set.seed(127)
    (data <- dodSim(tau=tau, d.prime=d, sample.size=100, method="user"))
    expect_false(givesWarnings(dod(data[1, ], data[2, ])))
})

#################################
## d.prime0 misspecification:
test_that("arguments are specified correctly", {
    dprime <- 1
    set.seed(127)
    data <- dodSim(d.prime=dprime, sample.size=100, method="equi.prob")

    expect_false(inherits(try(dod(data[1, ], data[2, ], d.prime0=1)), "try-error"))
    expect_error(dod(data[1, ], data[2, ], d.prime0=c(1, 2)))
    expect_error(dod(data[1, ], data[2, ], d.prime0=-1))
    expect_error(dod(data[1, ], data[2, ], d.prime0=TRUE))
    expect_error(dod(data[1, ], data[2, ], d.prime0="1"))

    expect_error(dod(data[1, ], data[2, ], d.prime0=0,
                     alternative="simil"),
                 "'alternative' has to be 'difference'")
    expect_error(dod(data[1, ], data[2, ], alternative="simil"),
                 "'alternative' has to be 'difference'")

    ## conf.level misspecification:
    expect_error(dod(data[1, ], data[2, ], conf.level=c(.9, .95)))
    expect_error(dod(data[1, ], data[2, ], conf.level=0))
    expect_error(dod(data[1, ], data[2, ], conf.level=1))

    ## statistic misspecification:
    expect_true(!inherits(try(dod(data[1, ], data[2, ], stat="like")), "try-error"))
    expect_true(!inherits(try(dod(data[1, ], data[2, ], stat="Pearson")), "try-error"))
    expect_true(!inherits(try(dod(data[1, ], data[2, ], stat="Wald")), "try-error"))
    expect_true(!inherits(try(dod(data[1, ], data[2, ], stat="Wilcox")), "try-error"))

    expect_error(dod(data[1, ], data[2, ], stat="FAIL"))

    expect_error(dod(data[1, ], data[2, ], stat="Wilcox", d.prime0=1),
                 "Wilcoxon statistic only available with d.prime0 = 0")

    ## 'alternative'-arg misspecification:
    expect_true(!inherits(try(dod(data[1, ], data[2, ],
                                  alternative="difference")), "try-error"))
    expect_true(!inherits(try(dod(data[1, ], data[2, ],
                                  alternative="greater")), "try-error"))
    expect_true(!inherits(try(dod(data[1, ], data[2, ], d.prime0=1,
                                  alternative="similarity")), "try-error"))
    expect_true(!inherits(try(dod(data[1, ], data[2, ], d.prime0=2,
                                  alternative="less")), "try-error"))
    expect_true(!inherits(try(dod(data[1, ], data[2, ], d.prime0=2,
                                  alternative="two.sided")), "try-error"))
    expect_true(inherits(try(dod(data[1, ], data[2, ], alternative="ox"),
                             silent=TRUE), "try-error"))

    expect_error(dod(data[1, ], data[2, ], d.prime0=0,
                     alternative="two.sided"),
                 "'alternative' has to be 'difference' or 'greater' if 'd.prime0' is 0")
    expect_error(dod(data[1, ], data[2, ],
                     alternative="two.sided"),
                 "'alternative' has to be 'difference' or 'greater' if 'd.prime0' is 0")
    expect_error(dod(data[1, ], data[2, ],
                     alternative="less"),
                 "'alternative' has to be 'difference' or 'greater' if 'd.prime0' is 0")
    expect_error(dod(data[1, ], data[2, ],
                     alternative="simil"),
                 "'alternative' has to be 'difference' or 'greater' if 'd.prime0' is 0")

    expect_error(dod(data[1, ], data[2, ], alternative="FAIL"))
})

test_that("Control settings work as intended", {
    ## Optimizer control settings:
    dprime <- 1
    set.seed(127)
    data <- dodSim(d.prime=dprime, sample.size=100, method="equi.prob")

    ctrl <- dodControl(optCtrl=list(trace=TRUE))
    expect_output(dod(data[1, ], data[2, ], control=ctrl),
                  "  0:     346.03849:  1.00000  2.00000  3.00000  1.00000")
    ## Faulty optimizer arguments:
    ctrl <- dodControl(optCtrl=list(verbose=TRUE))
    expect_warning(dod(data[1, ], data[2, ], control=ctrl),
                   "unrecognized control element")
    ## Force gradient warning from small gradient:
    ctrl <- dodControl(grad.tol=1e-10)
    expect_warning(dod(data[1, ], data[2, ], control=ctrl),
                   "Estimation failed")

    ## dodControl():
    expect_error(dodControl(grad.tol=-1e-10))
    expect_error(dodControl(grad.tol=Inf))
    expect_error(dodControl(grad.tol=NA))
    expect_error(dodControl(grad.tol=-Inf))

    ## Control settings:
    expect_true(inherits(dodControl(), "dodControl"))
    expect_error(dod(data[1, ], data[2, ], control=list()),
                 "Specify 'control' with dodControl()")
    expect_error(dod(data[1, ], data[2, ], control=NULL),
                 "Specify 'control' with dodControl()")
    expect_error(dod(data[1, ], data[2, ], control=NA),
                 "Specify 'control' with dodControl()")
    expect_true(!inherits(try(dod(data[1, ], data[2, ], control=dodControl())),
                          "try-error"))

    ## Turning off vcov and gradient computations:
    ctrl <- dodControl(get.vcov=FALSE)
    fm <- dod(data[1, ], data[2, ], control=ctrl)
    expect_true(is.na(coef(fm)[2]) && !is.null(fm$gradient))

    ctrl <- dodControl(get.grad=FALSE)
    fm <- dod(data[1, ], data[2, ], control=ctrl)
    expect_true(is.na(coef(fm)[2]) && is.null(fm$gradient))

    ## Not testing arguments:
    ctrl <- dodControl(test.args=FALSE)
    fm <- dod(data[1, ], data[2, ], d.prime0=-1, control=ctrl)
    expect_equal(fm$d.prime0, -1)
    data2 <- data
    data2[, 2] <- 0

    expect_error(
        givesWarnings(dod(data2[1, ], data2[2, ], control=ctrl,
                          stat="Wil"))
        )
})

test_that("warning works in get_tau()", {
    expect_warning(get_tau(1, ncat=4.000002),
                   "non-integer 'ncat': 4")
    expect_false(givesWarnings(get_tau(1, ncat=4.0000002)))
})

## same
## diff
## d.prime0
## conf.level
## statistic
## alternative
## control
## ... ## args to nlminb.

##################################################################
