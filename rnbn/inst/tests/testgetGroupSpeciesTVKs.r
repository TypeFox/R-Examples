context("Test getGroupSpeciesTVKs")
test_that("Returns a list of TVKS", {
    expect_true(length(getGroupSpeciesTVKs('reptile'))>10)    
})
