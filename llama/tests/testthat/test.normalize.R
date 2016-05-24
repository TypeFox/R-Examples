test_that("normalize normalizes", {
    normalized = normalize(data.frame(a = 1:10, b = 11:20))

    expect_equal(min(normalized$features$a), -1)
    expect_equal(max(normalized$features$a), 1)
    expect_equal(median(normalized$features$a), 0)
    expect_equal(min(normalized$features$b), -1)
    expect_equal(max(normalized$features$b), 1)
    expect_equal(median(normalized$features$b), 0)
})

test_that("normalize returns metadata", {
    normalized = normalize(data.frame(a = 1:10, b = 11:20))

    expect_equal(as.vector(normalized$meta$minValues), c(1, 11))
    expect_equal(as.vector(normalized$meta$maxValues), c(10, 20))
})

test_that("normalize takes metadata", {
    normalized = normalize(data.frame(a = 1:10, b = 11:20),
        list(minValues=data.frame(a=0, b=0), maxValues=data.frame(a=50, b=50)))

    expect_equal(min(normalized$features$a), -0.96)
    expect_equal(max(normalized$features$a), -0.6)
    expect_equal(median(normalized$features$a), -0.78)
    expect_equal(min(normalized$features$b), -0.56)
    expect_equal(max(normalized$features$b), -0.2)
    expect_equal(median(normalized$features$b), -0.38)
})
