context("UNFv6: Complex vectors")
test_that("Complex vectors optionally represented as A,iB", {
    expect_equal(unf6(0+1i, complex_as_character = FALSE)$unf, unf6("+0.e+,i+1.e+")$unf)
})
test_that("Complex vectors optionally treated as character", {
    expect_equal(unf6(0+1i, complex_as_character = TRUE)$unf, unf6("0+1i")$unf)
})
