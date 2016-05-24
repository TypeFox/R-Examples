require(tm)
test_that("Integration with tm", {
    expect_true(length(ngramrr("hello", char = TRUE, ngmax = 3, ngmin = 2)) == 7)
    expect_true(length(ngramrr("hello", char = TRUE, ngmax = 3)) == 12)
})
