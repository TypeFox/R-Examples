context("UNFv6: Matrices")
test_that("Matrix treated as data.frame", {
    expect_equal(unf(matrix(1:6, nrow=3), version = 6),
                 unf(as.data.frame(matrix(1:6, nrow=3)), version = 6))
})
test_that("Column order irrelevant", {
    expect_equal(unf(matrix(1:6, nrow=3), version = 6)$unf,
                 unf(matrix(1:6, nrow=3)[,2:1], version = 6)$unf)
})
test_that("Row order relevant", {
    expect_false(identical(unf(matrix(1:6, nrow=3), version = 6)$unf,
                           unf(matrix(1:6, nrow=3)[3:1,], version = 6)$unf))
})
test_that("Subsetting relevant", {
    expect_false(identical(unf(matrix(1:6, nrow=3), version = 6)$unf,
                           unf(matrix(1:6, nrow=3)[1:2,], version = 6)$unf))
})
