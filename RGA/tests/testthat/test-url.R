context("API request URL")

url <- get_url(path = c("management", "accounts"), query = list(start.index = 1, max.results = 1000))

test_that("URL class", {
    expect_is(url, "character")
})

test_that("URL length", {
    expect_equal(length(url), 1L)
})

test_that("URL match", {
    expect_match(url, "www.googleapis.com/analytics")
    expect_match(url, "management")
    expect_match(url, "accounts")
    expect_match(url, "start-index")
})
