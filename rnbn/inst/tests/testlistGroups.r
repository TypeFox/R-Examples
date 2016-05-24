context("Test listGroups")
groups <- listGroups()
test_that("List groups", {
    expect_true(is.data.frame(groups))
    expect_true(nrow(groups)>100)
})
