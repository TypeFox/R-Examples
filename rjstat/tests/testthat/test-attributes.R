context("Attributes")

dataset <- readLines("dataset.json")

test_that("attributes are correct", {
    expect_equal(attr(fromJSONstat(dataset)[[1]], "source"),
                 "Random data")
    expect_equal(attr(fromJSONstat(dataset)[[1]], "updated"),
                 "2014-09-29")
    expect_match(toJSONstat(fromJSONstat(dataset)),
                 "\"source\":\"Random data\"")
    expect_match(toJSONstat(fromJSONstat(dataset)),
                 "\"updated\":\"2014-09-29\"")
    expect_null(attr(fromJSONstat(dataset)[[2]], "source"))
    expect_null(attr(fromJSONstat(dataset)[[2]], "updated"))
})
