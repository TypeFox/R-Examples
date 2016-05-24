context("Test listVCs")
VCs <- listVCs()
test_that("List VCs", {
    expect_true(is.data.frame(VCs))
    expect_true(nrow(VCs)==112)
})
