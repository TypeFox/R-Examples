context("UNFv6: Missing Values")

test_that("Missing values calculated correctly", {
    expect_equal(unf6(NA)$unf, "cJ6AyISHokEeHuTfufIqhg==")
})

test_that("Nonfinites optionally treated as NA", {
    expect_equal(unf6(NaN, nonfinites_as_missing=TRUE)$unf, "cJ6AyISHokEeHuTfufIqhg==")
    expect_equal(unf6(Inf, nonfinites_as_missing=TRUE)$unf, "cJ6AyISHokEeHuTfufIqhg==")
    expect_equal(unf6(-Inf, nonfinites_as_missing=TRUE)$unf, "cJ6AyISHokEeHuTfufIqhg==")
})
