context("chart")

slug <- "obama-job-approval"
chartu <- pollstr_chart(slug, convert = FALSE)

test_that("API returns data", {
    expect_true(length(chartu) > 0)    
})

chart <- pollstr_chart(slug)

test_that("chart data is in the correct format", {
    expect_is(chart, "pollstr_chart")
    expect_equal(length(chart), 11)
})

