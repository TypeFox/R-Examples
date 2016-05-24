context("UNFv6: Characters")
test_that("Examples from original R package documentation", {
    expect_equal(unf6(c('test','1','2','3'))$unf, "fH4NJMYkaAJ16OWMEE+zpQ==")
})

test_that("Tails of long characters irrelevant", {
    lorem_ipsum <- c("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.", "Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.")
    expect_equal(unf6(lorem_ipsum), unf6(substring(lorem_ipsum, 1, 128)))
})

test_that("Tails of long characters optionally relevant", {
    lorem_ipsum <- c("Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.", "Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum.")
    expect_false(identical(unf6(lorem_ipsum), unf6(lorem_ipsum, characters=200, truncation=256)))
})

test_that("Numerics stored as character not same as numeric", {
    expect_false(unf6(c('1','2','3'))$unf == unf6(1:3)$unf)
})

test_that("Numerics stored as factors same as numeric stored as character", {
    expect_equal(unf6(c('1','2','3'))$unf, unf6(factor(c('1','2','3')))$unf)
})

test_that("truncation less than characters throws error", {
    expect_error(unf6(c('1','2','3'), characters = 128, truncation=5))
})
