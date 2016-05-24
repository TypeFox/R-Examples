context("sparsebnData")

dat <- generate_fixed_data_frame()
dat_na <- generate_na_data_frame()

test_that("sparsebnData constructor fails if input not a list, data.frame, or matrix", {
    expect_error(sparsebnData(1L))
    expect_error(sparsebnData(pi))
    expect_error(sparsebnData(rep(1,5)))
})

test_that("sparsebnData constructor fails if type not specified or improperly specified", {
    expect_error(sparsebnData(x = dat), regexp = "type")
    expect_error(sparsebnData(x = dat, type = "contns"), regexp = "Invalid \'type\'")
    expect_error(sparsebnData(x = dat, type = "dicsrt"), regexp = "Invalid \'type\'")

    ### Check that partial matching works OK
    expect_error(sparsebnData(x = dat, type = "c"), NA)
    expect_error(sparsebnData(x = dat, type = "cont"), NA)
    expect_error(sparsebnData(x = dat, type = "d"), NA)
    expect_error(sparsebnData(x = dat, type = "disc"), NA)
})

test_that("sparsebnData constructor issues a warning if data has missing values, but does not fail", {
    expect_warning(sparsebnData(x = dat_na, type = "continuous"), regexp = "Data contains [0-9]+")
})

test_that("print.sparsebnData functions properly", {
    expect_output(print(sparsebnData(x = dat, type = "continuous")), regexp = "5 total rows")
    expect_output(print(sparsebnData(x = dat, type = "continuous")), regexp = "Observational")
})
