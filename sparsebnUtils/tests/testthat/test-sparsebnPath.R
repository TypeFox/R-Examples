context("sparsebnPath")

test_that("sparsebnPath constructor fails if input not a list", {
    expect_error(sparsebnPath(1L))
    expect_error(sparsebnPath(pi))
    expect_error(sparsebnPath(rep(1,5)))
})

test_that("sparsebnPath constructor fails if components are not sparsebnFit", {
    expect_error(sparsebnPath(list(1L, 1L, 1L)))
    expect_error(sparsebnPath(list(pi, pi, pi)))
    expect_error(sparsebnPath(list(1:3, 4:6, 7:9)))
})
