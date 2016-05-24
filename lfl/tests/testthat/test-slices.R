test_that("slices 1", {
    from <- 34.5
    to <- 99.2
    n <- 1000
    res <- slices(from, to, n)
    expect_equal(min(res), from)
    expect_equal(max(res), to)
    expect_equal(length(res), n)
    expect_equal(sort(res), res)
})
