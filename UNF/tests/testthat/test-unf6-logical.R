context("UNFv6: Logicals")
test_that("Logical values treated as numeric 0/1s", {
    expect_equal(unf6(1), unf6(TRUE), label = "TRUE treated as numeric 1")
    expect_equal(unf6(0), unf6(FALSE), label = "FALSE treated as numeric 0")
})
