context("S4 class definition of \"sts\" and its extensions")

test_that("\"sts\" prototype is a valid object",
          expect_true(validObject(new("sts"))))

mysts <- sts(1:10, frequency = 4, start = c(1959, 2))

test_that("conversion from \"ts\" to \"sts\" works as expected", {
    myts <- ts(1:10, frequency = 4, start = c(1959, 2))
    expect_identical(as(myts, "sts"), mysts)
    ## this failed in surveillance 1.11.0 due to a wrong "start" calculation
})

test_that("if missing(observed), initialize-method copies slots", {
    mysts_updated <- initialize(mysts, epoch = 2:11)
    expect_identical(mysts_updated@epoch, 2:11)
    mysts_updated@epoch <- mysts@epoch
    expect_identical(mysts_updated, mysts)
    ## construct stsBP from existing "sts" object
    mystsBP <- new("stsBP", mysts,
                   ci = array(NA_real_, c(10,1,2)),
                   lambda = array(NA_real_, c(10,1,1)))
    expect_identical(as(mystsBP, "sts"), mysts)
})

test_that("different initializations of \"stsBP\" work as expected", {
    mystsBP <- new("stsBP", observed = 1:10, freq = 4, start = c(1959, 2),
                   ci = array(NA_real_, c(10,1,2)),
                   lambda = array(NA_real_, c(10,1,0)))
    expect_identical(mystsBP, as(mysts, "stsBP"))
})

test_that("different initializations of \"stsNC\" work as expected", {
    mystsNC <- new("stsNC", observed = 1:10, freq = 4, start = c(1959, 2),
                   pi = array(NA_real_, c(10,1,2)),
                   SR = array(NA_real_, c(10,0,0)))
    expect_identical(mystsNC, as(mysts, "stsNC"))
})

test_that("sts(..., population) sets the populationFrac slot", {
    ## for sts() construction, "population" is an alias for "populationFrac" 
    ## (the internal slot name), introduced in the space-time JSS paper
    sts1 <- sts(cbind(1:3, 11:13), population = c(10, 20))
    sts2 <- sts(cbind(1:3, 11:13), populationFrac = c(10, 20))
    expect_identical(sts1, sts2)
})
