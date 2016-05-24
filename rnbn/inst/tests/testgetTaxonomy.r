context("Test getTaxonomy")
test_that("Get taxonomy", {
    expect_true(is.data.frame(getTaxonomy("NHMSYS0000528028")))
    expect_true(nrow(getTaxonomy("NHMSYS0000528028")) == 6)
})
