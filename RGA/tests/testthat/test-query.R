context("Reporting API request query")

test_that("Ommit empty fields", {
    query <- build_query(profile.id = 0, sort = "", filters = NA)
    expect_null(query$sort)
    expect_null(query$filters)
})

test_that("Dates convert to a character", {
    query <- build_query(profile.id = 0, end.date = Sys.Date())
    expect_equal(query$end.date, as.character(Sys.Date()))
})

test_that("Strip white spaces", {
    query <- build_query(profile.id = 0, filters = "ga:users > 1000")
    expect_false(grepl(query$filters, " "))
})

test_that("Collapse all fields", {
    query <- build_query(profile.id = 0, metrics = c("ga:users", "ga:sessions"), dimensions = c("ga:date", "ga:hour"))
    expect_equal(length(query$metrics), 1L)
    expect_equal(length(query$dimensions), 1L)
})
