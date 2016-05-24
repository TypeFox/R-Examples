context("UNFv6: Lists")
test_that("Variable order irrelevant", {
    expect_equal(unf(list(1:3,4:6,7:9), version = 6)$unf,
                 unf(list(7:9,1:3,4:6), version = 6)$unf,
                 "ukDZSJXck7fn4SlPJMPFTQ==")
})
test_that("Sort order relevant", {
    expect_false(identical(unf(list(1:3,4:6,7:9), version = 6),
                           unf(list(rev(1:3),4:6,7:9), version = 6)))
})
test_that("Subsetting relevant", {
    expect_false(identical(unf(list(1:3,4:6,7:9), version = 6),
                           unf(list(1:2,4:5,7:8), version = 6)))
})
test_that("Dataframes and lists equivalent", {
    expect_equal(unf(data.frame(1:3,4:6,7:9), version = 6)$unf,
                 unf(list(1:3,4:6,7:9), version = 6)$unf,
                 "ukDZSJXck7fn4SlPJMPFTQ==")
})
