context("Test runnbnurl")
test_that("runnbnurl gets a response", {
    expect_that(length(runnbnurl(service="obs", tvks="NBNSYS0000002987", datasets="GA000373", startYear="1990", endYear="2010")) > 0, is_true()) 
})

    