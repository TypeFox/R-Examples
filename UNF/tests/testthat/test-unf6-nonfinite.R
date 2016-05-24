context("UNFv6: Non-finites")
test_that("Nonfinites handled per specification", {
    expect_equal(unf6("+inf"), unf6(Inf))
    expect_equal(unf6("-inf"), unf6(-Inf))
    expect_equal(unf6("+nan"), unf6(NaN))
})
test_that("Nonfinites optionally treated as missing", {
    expect_equal(unf6(NA), unf6(Inf, nonfinites_as_missing=TRUE))
    expect_equal(unf6(NA), unf6(-Inf, nonfinites_as_missing=TRUE))
    expect_equal(unf6(NA), unf6(NaN, nonfinites_as_missing=TRUE))
})
