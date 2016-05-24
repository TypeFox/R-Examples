context("Values")

dataset <- readLines("dataset.json")

test_that("values are correct", {
    expect_equal(fromJSONstat(dataset)[[1]]$value,
                 c(1.23456789, 2.3456789, 3.456789, 4.56789))
    expect_equal(fromJSONstat(dataset)[[2]]$value, c(NA, 2, NA, 4))
    d <- data.frame(V1 = rev(letters), value = 1:26)
    expect_equal(fromJSONstat(toJSONstat(d))[[1]]$value, 26:1)
})
