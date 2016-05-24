context("Test listDatasets")
lists<-listDatasets()
test_that("List datasets", {
    expect_true(is.data.frame(lists))
    expect_true(nrow(lists)>700)
})
